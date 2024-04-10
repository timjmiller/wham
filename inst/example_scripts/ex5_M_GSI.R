# WHAM example 5: Ecov and age-year effects on natural mortality

is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo
#by default do not perform bias-correction
if(!exists("basic_info")) basic_info <- NULL

library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)

# create directory for analysis, e.g.
if(!exists("write.dir")) write.dir <- tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)

# specify models:
# Model M_model       M_re  Ecov  Ecov-M link
# m1    ---           ---   ar1   ---
# m2    ---           ---   ar1   linear
# m3    ---           ---   ar1   quadratic
# m4    ---           ar1_a ar1   ---
# m5    ---           ar1_y ar1   ---
# m6    ---           2dar1 ar1   ---
# m7    age-specific  ---   ar1   ---
# m8    waa           ---   ar1   ---
# m9    constant      ---   ar1   ---
# m10   constant      ar1_a ar1   ---
# m11   constant      ar1_y ar1   ---
# m12   constant      2dar1 ar1   ---
# m13   constant      2dar1 ar1   linear
# m14   constant      2dar1 ar1   quadratic
Ecov_how <- paste0(
  c("none",rep("",2), rep("none", 9), rep("", 2)),
  c("", rep("lag-0-",2), rep("",9), rep("lag-0-",2)),
  c("", "linear", "poly-2", rep("",9), "linear", "poly-2"))

mean_model <- c(rep("fixed-M",6), "estimate-M", "weight-at-age", rep("estimate-M",6))
age_specific <- c(rep(NA,6),TRUE, NA, rep(FALSE, 6))

df.mods <- data.frame(M_model = c(rep("---",6),"age-specific","weight-at-age",rep("constant",6)),
                      mean_model = mean_model,
                      age_specific = age_specific,
                      M_re = c(rep("none",3),"ar1_a","ar1_y","ar1_ay",rep("none",3),"ar1_a", "ar1_y",rep("ar1_ay",3)),
                      Ecov_process = rep("ar1",14),
                      Ecov_how = Ecov_how, stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

mods <- vector("list",n.mods)
mods_proj <- vector("list",n.mods)
for(m in 1:n.mods){
  # set up environmental covariate data and model options
  # see ?prepare_wham_input
  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    M_how = array(df.mods$Ecov_how[m],c(1,1,asap3[[1]]$dat$n_ages,1))) # n_Ecov x n_stocks x n_ages x n_regions

  mean_map <- NULL
  if(df.mods$mean_model[m] == "estimate-M"){
    if(df.mods$age_specific[m]) mean_map[1,1,] <- 1:asap3[[1]]$dat$n_ages
    else mean_map[1,1,] <- 1
  }
  M <- list(
    mean_model = df.mods$mean_model[m],
    re_model = matrix(df.mods$M_re[m], 1,1),
    means_map = mean_map
  )
  if(df.mods$mean_model[m] == "estimate-M" & !df.mods$age_specific[m]) M$initial_means = array(0.28, c(1,1,asap3[[1]]$dat$n_ages)) #n_stocks x n_regions x n_ages

  # Generate wham input from ASAP3 and Ecov data
  print(paste("m:",m))
  input <- prepare_wham_input(asap3, recruit_model = 2,
    model_name = paste0("m",m,": ", df.mods$mean_model[m]," + GSI link: ",df.mods$Ecov_how[m]," + M RE: ", df.mods$M_re[m]),
    ecov = ecov,
    selectivity=list(model=rep("logistic",6),
      initial_pars=c(rep(list(c(3,3)),4), list(c(1.5,0.1), c(1.5,0.1))),
      fix_pars=c(rep(list(NULL),4), list(1:2, 1:2))),
    NAA_re = list(sigma='rec+1',cor='iid'),
    M=M,
    age_comp = "logistic-normal-pool0", basic_info = basic_info)

  # Fit model
  mods[[m]] <- fit_wham(input, do.retro=T, do.osa=F) # no OSA residuals to save time

  # Save model
  saveRDS(mods[[m]], file=paste0(df.mods$Model[m],".rds"))

  # If desired, plot output in new subfolder
  # plot_wham_output(mod=mod, out.type='html')

  # If desired, do projections
  # mod_proj[[i]] <- project_wham(mod)
  # saveRDS(mod_proj[[m]], file=paste0(df.mods$Model[m],"_proj.rds"))
}

# get convergence info
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)

# # make labeling prettier
# df.mods$Ecov_how[df.mods$Ecov_how=="none"] <- "---"
# df.mods$M_re[df.mods$M_re=="none"] = "---"
# colnames(df.mods)[2] = "M_est"

# only get AIC and Mohn's rho for converged models
df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
is_conv <- df.mods$conv & df.mods$pdHess
which(is_conv) # 1, 2, 5, 8, 9, 11, 12
mods2 <- mods[is_conv]
#mods2[not_conv] <- NULL
df.aic.tmp <- as.data.frame(compare_wham_models(mods2, table.opts=list(sort=FALSE, calc.rho=TRUE))$tab)
df.aic <- df.aic.tmp[FALSE,]
ct = 1
for(i in 1:n.mods){
  if(!is_conv[i]){
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
for(m in which(is_conv)){
  plot_wham_output(mod=mods[[m]], dir.main=file.path(write.dir,paste0("m",m)))
}

# save results table
write.csv(df.mods, file="ex5_table.csv",quote=F, row.names=F)

# compare models which(is_conv) 1, 2, 5, 8, 9, 11, 12
compare_wham_models(mods[which(is_conv)], do.table=FALSE, plot.opts=list(return.ggplot=F))
compare_wham_models(mods[which(is_conv)], do.table=FALSE, plot.opts=list(return.ggplot=F, which=6, M.age=5))
compare_wham_models(mods[which(is_conv)], do.table=FALSE, plot.opts=list(return.ggplot=F, which=6, M.age=4))

# ---------------------------------------------------------
# plot all MAA together in one giant 16-panel plot
years = mods[[1]]$years
n_years = length(years)
n_ages = mods[[1]]$env$data$n_ages
ages <- 1:n_ages

ecov_link <- df.mods$Ecov_how
ecov_link[ecov_link=="none"] = "no"
M_mod <- sapply(mods, function(x) x$input$options$M$mean_model)
M_re <- sapply(mods, function(x) x$input$options$M$re_model)
M_re[M_re=="none"] = "no"
df.MAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+3))
colnames(df.MAA) <- c(paste0("Age_",1:n_ages),"Year","Model","pdHess")
for(i in 1:n.mods){
  tmp = as.data.frame(mods[[i]]$rep$MAA[1,1,,])
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

png(filename = file.path(write.dir, "MAA.png"), width = 15, height = 15, res = 100, units='in')
    print(ggplot(df.plot2, aes(x=Year, y=Age)) +
      geom_tile(aes(fill=logM, alpha=factor(pdHess))) +
      scale_alpha_discrete(range=c(0.4,1), guide="none") +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw() +
      facet_wrap(~Model, nrow=4, dir="v") +
      scale_fill_viridis() #+ 
      # theme(strip.text = element_text(
      # size = 6)
    )
dev.off()

# add projections
mods_proj <- vector("list",n.mods)
tofit <- which(is_conv)
for(m in tofit){
  mods_proj[[m]] <- project_wham(mods[[m]], MakeADFun.silent=TRUE)
}

# make sure NLL doesn't change with projections
diff_nll <- sapply(tofit, function(x) mods[[m]]$opt$obj - mods_proj[[m]]$fn()) 
diff_nll #~0

