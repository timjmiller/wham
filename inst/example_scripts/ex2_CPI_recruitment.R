is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo
library(dplyr)
#by default do not perform bias-correction
if(!exists("basic_info")) basic_info <- NULL

# create directory for analysis, E.g.,
#write.dir <- "/path/to/save/output"
if(!exists("write.dir")) write.dir <- tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples, "CPI.csv"), header=T)

Ecov_how <- paste0(
  c("none", "controlling-", "none", "limiting-", "limiting-", "controlling-", "controlling-"), 
  c("", "lag-1-", "", rep("lag-1-",4)),
  c("", "linear", "", rep("linear", 4)))

df.mods <- data.frame(Recruitment = c(2,2,3,3,3,3,4),
                      Ecov_process = c(rep("rw",4),rep("ar1",3)),
                      Ecov_how = Ecov_how, stringsAsFactors=FALSE)

# specify 7 models:
# Model  Recruit_mod  Ecov_mod     Ecov_how
#    m1       Random        rw          ---
#    m2       Random        rw  Controlling
#    m3     Bev-Holt        rw          ---
#    m4     Bev-Holt        rw     Limiting
#    m5     Bev-Holt       ar1     Limiting
#    m6     Bev-Holt       ar1  Controlling
#    m7       Ricker       ar1  Controlling

n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

mods <- list()
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
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    recruitment_how = matrix(df.mods$Ecov_how[m],1,1)) #controlling  (dens-indep mortality), limiting (carrying capacity), lethal (threshold), masking (metabolism/growth), directive (behavior)

  # (not used in this vignette) can set ecov = NULL to fit model without Ecov data
  if(is.na(df.mods$Ecov_process[m])) ecov = NULL

  # Generate wham input from ASAP3 and Ecov data
  input <- prepare_wham_input(asap3, recruit_model = df.mods$Recruitment[m],
                              model_name = "Ex 2: SNEMA Yellowtail Flounder with CPI effects on R",
                              ecov = ecov,
                              NAA_re = list(sigma="rec+1", cor="iid"),
                              age_comp = "logistic-normal-pool0", basic_info = basic_info) # logistic normal pool 0 obs

  # Selectivity = logistic, not age-specific as in ex1
  #   2 pars per block instead of n.ages
  #   sel pars of indices 4/5 fixed at 1.5, 0.1 (specified via neg phase in ex2_SNEMAYT.dat)
  input$par$logit_selpars[1:4,7:8] <- 0 # last 2 rows will not be estimated (mapped to NA)

  # Fit model
  # collect fit models into a list
  mods[[m]] <- fit_wham(input, do.retro=TRUE, do.osa=TRUE)

  # Save model
  saveRDS(mods[[m]], file=file.path(write.dir, paste0(df.mods$Model[m],".rds")))

  # Plot output in new subfolder
  plot_wham_output(mod=mods[[m]], dir.main=file.path(write.dir,df.mods$Model[m]))
}


# check convergence of all models
vign2_conv <- lapply(mods, function(x) capture.output(check_convergence(x)))
for(m in 1:n.mods) cat(paste0("Model ",m,":"), vign2_conv[[m]], "", sep='\n')

# make results table prettier
df.mods$Recruitment <- dplyr::recode(df.mods$Recruitment, `2`='Random', `3`='Bev-Holt', `4`='Ricker')
df.mods$Ecov_how <- c("---", "Controlling", "---", "Limiting", "Limiting", "Controlling", "Controlling")

# get convergence info
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
df.mods$Rsig <- sapply(mods, function(x) round(exp(x$parList$log_NAA_sigma[1]),3))

# only get AIC and Mohn's rho for converged models
not_conv <- !df.mods$conv | !df.mods$pdHess
mods2 <- mods
mods2[not_conv] <- NULL
res <- compare_wham_models(mods2, table.opts=list(sort=FALSE, calc.rho=T))
df.aic.tmp <- as.data.frame(res$tab)
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

