# WHAM example 6: Numbers-at-age options

# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref='naa')
# devtools::load_all()
library(dplyr)
library(wham)

# create directory for analysis, e.g.
# write.dir <- "/path/to/save/ex2" on linux/mac
write.dir <- "/home/bstock/Documents/wham/sandbox/ex6_NAA"
if(!exists("write.dir")) write.dir = getwd()
dir.create(write.dir)
setwd(write.dir)

wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex2_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)

# confirm you are in the working directory and it has the ex2_SNEMAYT.dat file
list.files()

asap3 <- read_asap3_dat("ex2_SNEMAYT.dat")

# specify models:
# Model NAA_cor NAA_sigma
#    m1     ---       ---
#    m2     iid       rec
#    m3   ar1_a       rec
#    m4   ar1_y       rec
#    m5   2dar1       rec
#    m6     iid     rec+1
#    m7   ar1_a     rec+1
#    m8   ar1_y     rec+1
#    m9   2dar1     rec+1
df.mods <- data.frame(NAA_cor = c('---','iid','ar1_a','ar1_y','2dar1','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('---',rep("rec",4),rep("rec+1",4)), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

for(m in 1:n.mods){
  NAA_list <- list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"])
  if(NAA_list$sigma == '---') NAA_list = NULL
  input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = "Ex 6: Numbers-at-age",
                              selectivity=list(model=rep("logistic",6),
                                               initial_pars=c(rep(list(c(3,3)),4), list(c(1.5,0.1), c(1.5,0.1))),
                                               fix_pars=c(rep(list(NULL),4), list(1:2, 1:2))),
                              NAA_re = NAA_list)

  # overwrite age comp model (all models use logistic normal)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # Full state-space model, abundance is the state vector
  input$data$use_NAA_re = 1
  input$data$random_recruitment = 0
  input$map = input$map[!(names(input$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))] # remove from map (i.e. estimate)
  input$map$log_R = factor(rep(NA, length(input$par$log_R)))
  input$random = c(input$random, "log_NAA")

  # Fit model
  btime = Sys.time()
  mod <- fit_wham(input, do.retro=T, do.osa=F)
  mod$runtime = round(difftime(Sys.time(), btime, units = "mins"),1)

  # Save model
  saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))

  # Plot output in new subfolder
  # plot_wham_output(mod=mod, out.type='html')

  # Do projections
  # mod_proj <- project_wham(mod)
  # saveRDS(mod_proj, file=paste0(df.mods$Model[m],"_proj.rds"))
}

# collect fit models into a list
mod.list <- paste0(df.mods$Model,".rds")
mods <- lapply(mod.list, readRDS)
# mods_proj <- lapply(paste0(df.mods$Model,"_proj.rds"), readRDS)
df.mods$na_sdrep <- sapply(mods, function(x) x$na_sdrep)
# df.mods$pdHess <- sapply(mods, function(x) x$sdrep$pdHess)
df.mods$Ecov_link <- c("---","linear","poly-2")[df.mods$Ecov_link+1]
df.mods$M_re[df.mods$M_re=="none"] = "---"
colnames(df.mods)[2] = "M_est"

# deal with m11 and m13 having NaN NLL or NA sdrep
df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
theNA <- which(is.na(df.mods$NLL) | is.na(df.mods$na_sdrep))
mods2 <- mods
mods2[theNA] <- NULL
df.aic.tmp <- as.data.frame(compare_wham_models(mods2, sort=FALSE, calc.rho=T)$tab)
df.aic <- df.aic.tmp[FALSE,]
ct = 1
for(i in 1:n.mods){
  if(i %in% theNA){
    df.aic[i,] <- rep(NA,5)
  } else {
    df.aic[i,] <- df.aic.tmp[ct,]
    ct <- ct + 1
  }
}
df.aic$AIC[df.mods$na_sdrep==TRUE | is.na(df.mods$na_sdrep)] <- NA
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL

# look at results table
df.mods

# # plot output for all models
# for(m in 1:n.mods){
#   plot_wham_output(mod=mods[[m]], dir.main=file.path(getwd(),paste0("m",m)), out.type='html')
# }

# save results table
write.csv(df.mods, file="ex5_table.csv",quote=F, row.names=F)

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
df.MAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
colnames(df.MAA) <- c(paste0("Age_",1:n_ages),"Year","Model")
for(i in 1:n.mods){
  tmp = as.data.frame(mods[[i]]$rep$MAA)
  tmp$Year <- years
  colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
  tmp$Model = paste0("m",i,": ", M_mod[i]," + ",ecov_link[i]," GSI link + ",M_re[i]," devs")
  df.MAA <- rbind(df.MAA, tmp)
}
df.plot <- df.MAA %>% tidyr::pivot_longer(-c(Year,Model),
          names_to = "Age",
          names_prefix = "Age_",
          names_ptypes = list(Age = integer()),
          values_to = "M")
df.plot2 <- dplyr::filter(df.plot, ! Model %in% c("m13: age-specific + poly-2 GSI link + 2D AR1 devs","m11: constant + linear GSI link + 2D AR1 devs"))
df.plot2$Model <- factor(as.character(df.plot2$Model), levels=unique(df.plot2$Model))
df.plot2$logM <- log(df.plot2$M)
df.plot2$logM[df.plot2$logM < -4] <- -4

png(filename = file.path(getwd(), paste0("MAA.png")), width = 10, height = 8, res = 100, units='in')
    print(ggplot2::ggplot(df.plot2, ggplot2::aes(x=Year, y=Age, fill=logM)) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_continuous(expand=c(0,0)) +
      ggplot2::scale_y_continuous(expand=c(0,0)) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~Model, nrow=5, dir="v") +
      viridis::scale_fill_viridis())
dev.off()

# # make sure NLL doesn't change with projections
# samenll <- mapply(function(x,y) ifelse(all.equal(x$opt$objective, y$opt$objective)==TRUE, TRUE, FALSE), mods, mods_proj)
# samenll
