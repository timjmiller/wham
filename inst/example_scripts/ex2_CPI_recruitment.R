# WHAM example 2: Cold Pool Index effect on SNEMA Yellowtail Flounder recruitment
# Replicate Miller et al 2016 results
#   adds environmental covariate (CPI, treated as rw)
#   uses 5 indices (ex 1: only 2 indices)
#   only fit to 1973-2011 data (ex 1: 1973-2016)
#   age compositions = 5, logistic normal pool zero obs (ex 1: 7, logistic normal missing zero obs)
#   selectivity = logistic (ex 1: age-specific)

# load wham
library(wham)
library(dplyr)

# create directory for analysis, e.g.
# write.dir <- "/path/to/save/ex2" on linux/mac
if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

# copy data files to working directory
wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex2_SNEMAYT.dat"), to=write.dir, overwrite=TRUE)
file.copy(from=file.path(wham.dir,"extdata","CPI.csv"), to=write.dir, overwrite=TRUE)

# confirm you are in the working directory and it has the ex2_SNEMAYT.dat and CPI.csv files
list.files()

# read asap3 data file and convert to input list for wham
asap3 <- read_asap3_dat("ex2_SNEMAYT.dat")

# load env covariate, CPI (Cold Pool Index)
env.dat <- read.csv("CPI.csv", header=T)

# specify 7 models:
# Model  Recruit_mod  Ecov_mod     Ecov_how
#    m1       Random        rw          ---
#    m2       Random        rw  Controlling
#    m3     Bev-Holt        rw          ---
#    m4     Bev-Holt        rw     Limiting
#    m5     Bev-Holt       ar1     Limiting
#    m6     Bev-Holt       ar1  Controlling
#    m7       Ricker       ar1  Controlling
df.mods <- data.frame(Recruitment = c(2,2,3,3,3,3,4),
                      Ecov_process = c(rep("rw",4),rep("ar1",3)),
                      Ecov_how = c(0,1,0,2,2,1,1), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

# run the 7 models
for(m in 1:n.mods){
  # set up environmental covariate data and model options
  # see ?prepare_wham_input
  ecov <- list(
    label = "CPI",
    mean = as.matrix(env.dat$CPI),
    logsigma = as.matrix(log(env.dat$CPI_sigma)),
    year = env.dat$Year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    lag = 1, # CPI in year t affects recruitment in year t+1
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    where = c("none","recruit")[as.logical(df.mods$Ecov_how[m])+1], # 'recruit' = CPI affects recruitment, 'none' = no effect (how = 0)
    how = df.mods$Ecov_how[m]) # 0 = no effect (but still fit Ecov to compare AIC), 1 = controlling (dens-indep mortality), 2 = limiting (carrying capacity), 3 = lethal (threshold), 4 = masking (metabolism/growth), 5 = directive (behavior)

  # (not used in this vignette) can set ecov = NULL to fit model without Ecov data
  if(is.na(df.mods$Ecov_process[m])) ecov = NULL

  # Generate wham input from ASAP3 and Ecov data
  input <- prepare_wham_input(asap3, recruit_model = df.mods$Recruitment[m],
                              model_name = "Ex 2: SNEMA Yellowtail Flounder with CPI effects on R",
                              ecov = ecov,
                              NAA_re = list(sigma="rec+1", cor="iid"),
                              age_comp = "logistic-normal-pool0") # logistic normal pool 0 obs

  # Selectivity = logistic, not age-specific as in ex1
  #   2 pars per block instead of n.ages
  #   sel pars of indices 4/5 fixed at 1.5, 0.1 (specified via neg phase in ex2_SNEMAYT.dat)
  input$par$logit_selpars[1:4,7:8] <- 0 # last 2 rows will not be estimated (mapped to NA)

  # Fit model
  mod <- fit_wham(input, do.retro=TRUE, do.osa=TRUE)

  # Save model
  saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))

  # Plot output in new subfolder
  plot_wham_output(mod=mod, dir.main=file.path(getwd(),df.mods$Model[m]), out.type='png')
}

# collect fit models into a list
mod.list <- paste0(df.mods$Model,".rds")
mods <- lapply(mod.list, readRDS)

# check convergence of all models
vign2_conv <- lapply(mods, function(x) capture.output(check_convergence(x)))
for(m in 1:n.mods) cat(paste0("Model ",m,":"), vign2_conv[[m]], "", sep='\n')

# make results table prettier
df.mods$Recruitment <- dplyr::recode(df.mods$Recruitment, `2`='Random', `3`='Bev-Holt', `4`='Ricker')
df.mods$Ecov_how <- dplyr::recode(df.mods$Ecov_how, `0`='---',`1`='Controlling', `2`='Limiting', `4`='Masking')

# get convergence info
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))

# only get AIC and Mohn's rho for converged models
not_conv <- !df.mods$conv | !df.mods$pdHess
mods2 <- mods
mods2[not_conv] <- NULL
df.aic.tmp <- as.data.frame(compare_wham_models(mods2, table.opts=list(sort=FALSE, calc.rho=T))$tab)
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
df.mods <- cbind(df.mods, df.aic)
df.mods <- df.mods[order(df.mods$dAIC, na.last=TRUE),]
df.mods[is.na(df.mods$AIC), c('dAIC','AIC','rho_R','rho_SSB','rho_Fbar')] <- "---"
rownames(df.mods) <- NULL

# look at results table
df.mods

# save results table
save("df.mods", file="vign2_res.RData")
write.csv(df.mods, file="vign2_res.csv",row.names=F,quote=F)
