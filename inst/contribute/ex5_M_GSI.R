library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
setwd(write.dir)
file.copy(from=file.path(pkg.dir,"inst","extdata","ex2_SNEMAYT.dat"), to=write.dir, overwrite=TRUE)
file.copy(from=file.path(pkg.dir,"inst","extdata","GSI.csv"), to=write.dir, overwrite=TRUE)
asap3 <- read_asap3_dat("ex2_SNEMAYT.dat")
env.dat <- read.csv("GSI.csv", header=T)

# specify models:
# Model M_model       M_re  Ecov  Ecov-M link
# m1    ---           ---   ar1   ---
# m2    ---           ---   ar1   linear
# m3    ---           ---   ar1   quadratic
# m4    age-specific  ---   ar1   ---
# m5    waa           ---   ar1   ---
# m6    constant      ---   ar1   ---
# m7    constant      ar1_y ar1   ---
# m8    constant      2dar1 ar1   ---
# m9    constant      ---   ar1   linear
# m10   constant      ---   ar1   quadratic
# m11   ---           2dar1 ar1   ---
# df.mods <- data.frame(M_model = c(rep("---",3),"age-specific","weight-at-age",rep("constant",6),"age-specific","age-specific",rep("constant",3),"---"),
#                       M_re = c(rep("none",6),"ar1_y","2dar1","none","none","2dar1","none","2dar1",rep("ar1_a",3),"2dar1"),
#                       Ecov_process = rep("ar1",17),
#                       Ecov_link = c(0,1,2,rep(0,5),1,2,1,2,2,0,1,2,0), stringsAsFactors=FALSE)
df.mods <- data.frame(M_model = c(rep("---",3),"age-specific","weight-at-age",rep("constant",5),"---"),
                      M_re = c(rep("none",6),"ar1_y","2dar1","none","none","2dar1"),
                      Ecov_process = rep("ar1",11),
                      Ecov_link = c(0,1,2,rep(0,5),1,2,0), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

for(m in 1:n.mods){
  # set up environmental covariate data and model options
  # see ?prepare_wham_input
  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    lag = 0, # GSI in year t affects M in same year
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    where = "M", # GSI affects natural mortality
    how = ifelse(df.mods$Ecov_link[m]==0,0,1), # 0 = no effect (but still fit Ecov to compare AIC), 1 = mean
    link_model = c(NA,"linear","poly-2")[df.mods$Ecov_link[m]+1])

  m_model <- df.mods$M_model[m]
  if(df.mods$M_model[m] == '---') m_model = "age-specific"
  if(df.mods$M_model[m] %in% c("constant","weight-at-age")) est_ages = 1
  if(df.mods$M_model[m] == "age-specific") est_ages = 1:asap3$dat$n_ages
  if(df.mods$M_model[m] == '---') est_ages = NULL
  M <- list(
    model = m_model,
    re = df.mods$M_re[m],
    est_ages = est_ages
  )
  if(m_model == "constant") M$initial_means = 0.28

  # Generate wham input from ASAP3 and Ecov data
  # input <- prepare_wham_input(asap3, recruit_model = 2,
  #                             model_name = "Ex 5: Yellowtail Flounder with GSI effects on M",
  #                             ecov = ecov,
  #                             M = M)
  input <- prepare_wham_input(asap3, recruit_model = 2,
                              model_name = paste0("m",m,": ", df.mods$M_model[m]," + ",c("no","linear","poly-2")[df.mods$Ecov_link[m]+1]," GSI link + ",df.mods$M_re[m]," devs"),
                              ecov = ecov,
                              selectivity=list(model=rep("logistic",6),
                                               initial_pars=c(rep(list(c(3,3)),4), list(c(1.5,0.1), c(1.5,0.1))),
                                               fix_pars=c(rep(list(NULL),4), list(1:2, 1:2))),
                              NAA_re = list(sigma='rec+1',cor='iid'),
                              M=M,
                              age_comp = "logistic-normal-pool0") # logistic normal pool 0 obs

  # Fit model
  mod <- fit_wham(input, do.retro=T, do.osa=F) # no OSA residuals to save time

  # Save model
  saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))

  # If desired, plot output in new subfolder
  # plot_wham_output(mod=mod, out.type='html')

  # If desired, do projections
  # mod_proj <- project_wham(mod)
  # saveRDS(mod_proj, file=paste0(df.mods$Model[m],"_proj.rds"))
}

# collect fit models into a list
mod.list <- paste0(df.mods$Model,".rds")
mods <- lapply(mod.list, readRDS)
# mods_proj <- lapply(paste0(df.mods$Model,"_proj.rds"), readRDS)

# get convergence info
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)

# make labeling prettier
df.mods$Ecov_link <- c("---","linear","poly-2")[df.mods$Ecov_link+1]
df.mods$M_re[df.mods$M_re=="none"] = "---"
colnames(df.mods)[2] = "M_est"

# only get AIC and Mohn's rho for converged models
df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
not_conv <- !df.mods$conv | !df.mods$pdHess
mods2 <- mods
mods2[not_conv] <- NULL
df.aic.tmp <- as.data.frame(compare_wham_models(mods2, table.opts=list(sort=FALSE, calc.rho=TRUE))$tab)
df.aic <- df.aic.tmp[FALSE,]
ct = 1
for(i in 1:n.mods){
  if(not_conv[i]){
    df.aic[i,] <- rep(NA,5)
  } else {
    df.aic[i,] <- df.aic.tmp[ct,]
    ct <- ct + 1
  }
}
df.aic[,1:2] <- format(round(df.aic[,1:2], 1), nsmall=1)
df.aic[,3:5] <- format(round(df.aic[,3:5], 3), nsmall=3)
df.aic[grep("NA",df.aic$dAIC),] <- "---"
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL

# look at results table
df.mods

# plot output for all models that converged
for(m in which(!not_conv)){
  plot_wham_output(mod=mods[[m]], dir.main=file.path(getwd(),paste0("m",m)), out.type='png')
}

# save results table
write.csv(df.mods, file="ex5_table.csv",quote=F, row.names=F)

# compare models 1, 2, 6, 8, 11
compare_wham_models(mods[c(1,2,6,8,11)], do.table=FALSE, plot.opts=list(return.ggplot=F))
compare_wham_models(mods[c(1,2,6,8,11)], do.table=FALSE, plot.opts=list(return.ggplot=F, which=6, M.age=5))
compare_wham_models(mods[c(1,2,6,8,11)], do.table=FALSE, plot.opts=list(return.ggplot=F, which=6, M.age=4))

# ---------------------------------------------------------
# plot all MAA together in one giant 16-panel plot
years = mods[[1]]$years
n_years = length(years)
n_ages = mods[[1]]$env$data$n_ages
ages <- 1:n_ages

ecov_link <- df.mods$Ecov_link
ecov_link[ecov_link=="---"] = "no"
M_mod <- c("constant","age-specific","weight-at-age")[sapply(mods, function(x) x$env$data$M_model)]
M_mod[sapply(mods, function(x) unique(x$env$data$M_est)) == 0] = "fixed"
M_re <- c("no","IID","AR1_a","AR1_y","2D AR1")[sapply(mods, function(x) x$env$data$M_re_model)]
df.MAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+3))
colnames(df.MAA) <- c(paste0("Age_",1:n_ages),"Year","Model","pdHess")
for(i in 1:n.mods){
  tmp = as.data.frame(mods[[i]]$rep$MAA)
  tmp$Year <- years
  colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
  tmp$Model = paste0("m",i,": ", M_mod[i]," + ",ecov_link[i]," GSI link + ",M_re[i]," devs")
  tmp$pdHess <- df.mods$pdHess[i]
  df.MAA <- rbind(df.MAA, tmp)
}
df.plot <- df.MAA %>% tidyr::pivot_longer(-c(Year,Model,pdHess),
          names_to = "Age",
          names_prefix = "Age_",
          names_transform = list(Age = as.integer),
          values_to = "M")
# df.plot2 <- dplyr::filter(df.plot, ! Model %in% c("m13: age-specific + poly-2 GSI link + 2D AR1 devs","m11: constant + linear GSI link + 2D AR1 devs"))
df.plot2 <- df.plot
df.plot2$Model <- factor(as.character(df.plot2$Model), levels=unique(df.plot2$Model))
df.plot2$logM <- log(df.plot2$M)
df.plot2$logM[df.plot2$logM < -4] <- -4

png(filename = file.path(getwd(), paste0("MAA.png")), width = 8.5, height = 5, res = 100, units='in')
    print(ggplot(df.plot2, aes(x=Year, y=Age)) +
      geom_tile(aes(fill=logM, alpha=factor(pdHess))) +
      scale_alpha_discrete(range=c(0.4,1), guide=FALSE) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw() +
      facet_wrap(~Model, nrow=4, dir="v") +
      scale_fill_viridis())
dev.off()

# add projections
mods_proj <- vector("list",n.mods)
tofit <- c(1:3,5:11)
for(m in tofit){
  mods_proj[[m]] <- project_wham(mods[[m]], MakeADFun.silent=TRUE)
}

# make sure NLL doesn't change with projections
samenll <- mapply(function(x,y) ifelse(all.equal(x$opt$objective, y$opt$objective)==TRUE, TRUE, FALSE), mods[tofit], mods_proj[tofit])
samenll

