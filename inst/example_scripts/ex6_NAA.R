# WHAM example 6: Numbers-at-age options

# As in example 1:
#   stock = SNEMA yellowtail flounder
#   2 indices
#   fit to 1973-2016 data
#   age compositions = 7, logistic normal don't pool zero obs

# As in example 2:
#   environmental effect on recruitment

# As in example 4:
#   selectivity = age-specific (fix ages = list(4:5,4,2:4))

# As in example 5:
#   Gulf Stream Index (GSI)

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(ggplot2)
library(tidyr)
library(dplyr)

# create directory for analysis, e.g.
# write.dir <- "/path/to/save/ex2" on linux/mac
if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex1_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)
file.copy(from=file.path(wham.dir,"extdata","GSI.csv"), to=write.dir, overwrite=FALSE)

# confirm you are in the working directory and it has the ex2_SNEMAYT.dat file
list.files()

asap3 <- read_asap3_dat("ex1_SNEMAYT.dat")
env.dat <- read.csv("GSI.csv", header=T)

# specify models:
#    Model NAA_cor NAA_sigma GSI_how
# 1     m1     ---       ---       0
# 2     m2     iid       rec       0
# 3     m3   ar1_y       rec       0
# 4     m4     iid     rec+1       0
# 5     m5   ar1_a     rec+1       0
# 6     m6   ar1_y     rec+1       0
# 7     m7   2dar1     rec+1       0
# 8     m8     iid       rec       2
# 9     m9   ar1_y       rec       2
# 10   m10     iid     rec+1       2
# 11   m11   ar1_a     rec+1       2
# 12   m12   ar1_y     rec+1       2
# 13   m13   2dar1     rec+1       2
df.mods <- data.frame(NAA_cor = c('---','iid','ar1_y','iid','ar1_a','ar1_y','2dar1','iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('---',rep("rec",2),rep("rec+1",4),rep("rec",2),rep("rec+1",4)),
                      GSI_how = c(rep(0,7),rep(2,6)), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

for(m in 1:n.mods){
  NAA_list <- list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"])
  if(NAA_list$sigma == '---') NAA_list = NULL

  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    lag = 1, # GSI in year t affects Rec in year t + 1
    process_model = 'ar1', # "rw" or "ar1"
    where = c("none","recruit")[as.logical(df.mods$GSI_how[m])+1], # GSI affects recruitment
    how = df.mods$GSI_how[m], # 0 = no effect (but still fit Ecov to compare AIC), 2 = limiting
    link_model = "linear")

  input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = "Ex 6: Numbers-at-age",
                              selectivity=list(model=rep("age-specific",3), re=c("none","none","none"), 
                                initial_pars=list(c(0.1,0.5,0.5,1,1,1),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,0.5,1,1,1,1)), 
                                fix_pars=list(4:6,4,3:6)),
                              NAA_re = NAA_list,
                              ecov=ecov,
                              age_comp = "logistic-normal-miss0") # logistic normal, treat 0 obs as missing

  # Fit model
  mod <- fit_wham(input, do.retro=T, do.osa=F)

  # Save model
  saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))
}

# collect fit models into a list
mod.list <- paste0(df.mods$Model,".rds")
mods <- lapply(mod.list, readRDS)

# get convergence info
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)

# make labeling prettier
df.mods$GSI_how <- c("---","---","Limiting")[df.mods$GSI_how+1]

# only get AIC and Mohn's rho for converged models
df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
not_conv <- !df.mods$conv | !df.mods$pdHess
mods2 <- mods
mods2[not_conv] <- NULL
df.aic.tmp <- as.data.frame(compare_wham_models(mods2, sort=FALSE, calc.rho=T)$tab)
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
mods[[1]]$env$data$recruit_model = 2 # m1 didn't actually fit a Bev-Holt
for(m in which(!not_conv)){
  plot_wham_output(mod=mods[[m]], dir.main=file.path(getwd(),paste0("m",m)), out.type='html')
}

# save results table
write.csv(df.mods, file="ex6_table.csv",quote=F, row.names=F)

# ---------------------------------------------------------
# plot all NAA_re in one multipanel plot
years = mods[[1]]$years[-1]
n_years = length(years)
n_ages = mods[[1]]$env$data$n_ages
ages <- 1:n_ages

plot.mods <- which(!not_conv)[-1] # m1 doesn't have NAA devs bc no stock-recruit function to predict rec
NAA_mod <- c("FE","RE: Recruit","RE: all NAA")[sapply(mods[plot.mods], function(x) x$env$data$n_NAA_sigma+1)]
NAA_cor <- c("IID","AR1_a","AR1_y","2D AR1")[sapply(mods[plot.mods], function(x) 4-sum(which(x$parList$trans_NAA_rho == 0)))]
GSI_how <- c("no GSI-Recruitment link","GSI-Recruitment link (limiting)")[as.numeric(factor(df.mods$GSI_how[plot.mods]))]
NAA_lab <- paste(NAA_mod,NAA_cor,sep=" + ") # duplicate by GSI, for facet_grid
df.NAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+3))
colnames(df.NAA) <- c(paste0("Age_",1:n_ages),"Year","GSI_how","NAA_lab")
for(i in 1:length(plot.mods)){
  tmp = as.data.frame(mods[[plot.mods[i]]]$rep$NAA_devs)
  tmp$Year <- years
  colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
  tmp$GSI_how = GSI_how[i]
  tmp$NAA_lab = NAA_lab[i]
  df.NAA <- rbind(df.NAA, tmp)
}
df.plot <- df.NAA %>% tidyr::pivot_longer(-c(Year,GSI_how,NAA_lab),
          names_to = "Age",
          names_prefix = "Age_",
          names_transform = list(Age = as.integer),
          values_to = "NAA_re")
NAA_lab_levels <- c("RE: Recruit + IID","RE: Recruit + AR1_y","RE: all NAA + IID","RE: all NAA + AR1_a","RE: all NAA + AR1_y","RE: all NAA + 2D AR1")
GSI_lab_levels <- c("no GSI-Recruitment link","GSI-Recruitment link (limiting)")
df.plot$GSI_how <- factor(as.character(df.plot$GSI_how), levels=GSI_lab_levels)
df.plot$NAA_lab <- factor(as.character(df.plot$NAA_lab), levels=NAA_lab_levels)
df.plot$NAA_re[df.plot$NAA_lab %in% c("RE: Recruit + IID","RE: Recruit + AR1_y") & df.plot$Age > 1] = 0

png(filename = file.path(getwd(), paste0("NAA_devs.png")), width = 8, height = 8.5, res = 200, units='in')
    print(ggplot(df.plot, ggplot2::aes(x=Year, y=Age)) +
      geom_tile(aes(fill=NAA_re)) +
      geom_label(aes(x=Year, y=Age, label=lab), size=5, alpha=1, #fontface = "bold",
          data=data.frame(Year=1976.5, Age=5.8, lab=df.mods$Model[plot.mods], NAA_lab=factor(NAA_lab, levels=NAA_lab_levels), GSI_how=factor(GSI_how, levels=GSI_lab_levels))) +                
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw() +
      facet_grid(rows=vars(NAA_lab), cols=vars(GSI_how), drop=F) +
      scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white", high = scales::muted("red")))
dev.off()

