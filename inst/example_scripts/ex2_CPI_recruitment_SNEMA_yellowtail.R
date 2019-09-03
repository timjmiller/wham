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
dir.create(write.dir)
setwd(write.dir)

# copy data files to working directory
wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex2_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)
file.copy(from=file.path(wham.dir,"extdata","CPI.csv"), to=write.dir, overwrite=FALSE)

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
  Ecov <- list(
    label = "CPI",
    mean = as.matrix(env.dat$CPI),
    sigma = as.matrix(env.dat$CPI_sigma),
    year = env.dat$Year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    lag = 1, # CPI in year t affects recruitment in year t+1
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    where = "recruit", # CPI affects recruitment
    how = df.mods$Ecov_how[m]) # 0 = no effect (but still fit Ecov to compare AIC), 1 = controlling (dens-indep mortality), 2 = limiting (carrying capacity), 3 = lethal (threshold), 4 = masking (metabolism/growth), 5 = directive (behavior)

  # (not used in this vignette) can set Ecov = NULL to fit model without Ecov data
  if(is.na(df.mods$Ecov_process[m])) Ecov = NULL

  # Generate wham input from ASAP3 and Ecov data
  input <- prepare_wham_input(asap3, recruit_model = df.mods$Recruitment[m],
                              model_name = "Ex 2: SNEMA Yellowtail Flounder with CPI effects on R",
                              Ecov = Ecov)

  # Builds off model m4 in example 1:
  #   full state-space model, logistic normal age-compositions

  # Age comp model = logistic normal pool obs (not multinomial, the default)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # Selectivity = logistic, not age-specific as in ex1
  #   2 pars per block instead of n.ages
  #   sel pars of indices 4/5 fixed at 1.5, 0.1 (specified via neg phase in ex2_SNEMAYT.dat)
  input$par$logit_selpars[1:4,7:8] <- 0 # last 2 rows will not be estimated (mapped to NA)

  # Full state-space model, abundance is the state vector
  input$data$use_NAA_re = 1
  input$data$random_recruitment = 0
  input$map = input$map[!(names(input$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
  input$map$log_R = factor(rep(NA, length(input$par$log_R)))
  input$random = c(input$random, "log_NAA","Ecov_re")

  # ---------------------------------------------------------
  ## Fit model
  mod <- fit_wham(input, do.retro=TRUE, do.osa=TRUE)

  # Save model
  saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))

  # Plot output in new subfolder
  plot_wham_output(mod=mod, dir.main=file.path(getwd(),df.mods$Model[m]), out.type='html')
}

# collect fit models into a list
mod.list <- paste0(df.mods$Model,".rds")
mods <- lapply(mod.list, readRDS)

# check convergence of all models
vign2_conv <- lapply(mods, function(x) capture.output(check_convergence(x)))
for(m in 1:n.mods) cat(paste0("Model ",m,":"), vign2_conv[[m]], "", sep='\n')

# calculate AIC and Mohn's rho
df.aic <- compare_wham_models(mods, sort=FALSE)$tab
df.mods <- cbind(df.mods, df.aic)

# make results prettier
rownames(df.mods) <- NULL
df.mods$Recruitment <- dplyr::recode(df.mods$Recruitment, `2`='Random', `3`='Bev-Holt', `4`='Ricker')
df.mods$Ecov_how <- dplyr::recode(df.mods$Ecov_how, `0`='---',`1`='Controlling', `2`='Limiting', `4`='Masking')
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
df.mods <- df.mods[order(df.mods$dAIC),]

# look at results table
df.mods

# save results table
save("df.mods", file="vign2_res.RData")
