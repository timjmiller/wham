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
file.copy(from=file.path(wham.dir,"extdata","ASAP_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)
file.copy(from=file.path(wham.dir,"extdata","CPI.csv"), to=write.dir, overwrite=FALSE)

# confirm you are in the working directory and it has the ASAP_SNEMAYT.dat and CPI.csv files
list.files()

# read asap3 data file and convert to input list for wham
asap3 <- read_asap3_dat("ASAP_SNEMAYT.dat")

# load env covariate, CPI (Cold Pool Index)
env.dat <- read.csv("CPI.csv", header=T)

# specify 7 models:
# Model  Recruit_mod  Ecov_mod     Ecov_how
#    m1       Random       ar1          ---
#    m2       Random       ar1  Controlling
#    m3     Bev-Holt       ar1          ---
#    m4     Bev-Holt       ar1     Limiting
#    m5     Bev-Holt        rw     Limiting
#    m6     Bev-Holt       ar1  Controlling
#    m7       Ricker       ar1  Controlling
df.mods <- data.frame(Recruitment = c(2,2,3,3,3,3,4),
                      Ecov_process = c('ar1','ar1','ar1','ar1','rw','ar1','ar1'),
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
    sigma = as.matrix(env.dat$CPI_sigma),
    year = env.dat$Year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    lag = 1, # CPI in year t affects recruitment in year t+1
    process_model = df.mods$Ecov_process[m], # c("rw", "ar1")
    where = "recruit", # CPI affects population via recruitment
    how = df.mods$Ecov_how[m]) # 0 = no effect (but still fit Ecov to compare AIC), 1 = controlling (dens-indep mortality), 2 = limiting (carrying capacity), 3 = lethal (threshold), 4 = masking (metabolism/growth), 5 = directive (behavior)

  # (not used in this vignette) can set ecov = NULL to fit model without ecov data
  if(is.na(df.mods$Ecov_process[m])) ecov = NULL 

  # generate wham input from ASAP3 and ecov data
  input <- prepare_wham_input(asap3, recruit_model = df.mods$Recruitment[m],
                              model_name = "SNEMA Yellowtail Flounder with CPI effects on R",
                              ecov = ecov)

  # builds off model m1 in example 1:
  #   SCAA, but with random effects for recruitment and index observation error variances fixed  

  # make one or more selectivity blocks with age-specific parameters
  age.specific = 1:3 # 3 age-specific blocks
  not.age.specific = (1:input$data$n_selblocks)[-age.specific]
  input = set_age_sel0(input, age.specific)
  input$par$logit_selpars[not.age.specific,c(1:input$data$n_ages,input$data$n_ages + 3:6)] = Inf
  input$par$logit_selpars[1,5] = Inf
  input$par$logit_selpars[2,4] = Inf
  input$par$logit_selpars[3,2] = Inf

  # Now redefine the map argument for the selectivity parameters to estimate only selectivity parameters without initial values at lower and upper bounds.
  input$map$logit_selpars = matrix(input$map$logit_selpars, input$data$n_selblocks, input$data$n_ages + 6)
  input$map$logit_selpars[is.infinite(input$par$logit_selpars)] = NA
  input$map$logit_selpars[!is.infinite(input$par$logit_selpars)] = 1:sum(!is.infinite(input$par$logit_selpars))
  input$map$logit_selpars = factor(input$map$logit_selpars)

  # full state-space model, abundance is the state vector
  input$data$use_NAA_re = 1
  input$data$random_recruitment = 0
  input$map = input$map[!(names(input$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
  input$map$log_R = factor(rep(NA, length(input$par$log_R)))
  input$random = "log_NAA"

  # fit model
  mod <- fit_wham(input, do.retro=TRUE, do.osa=TRUE)

  # check convergence
  check_convergence(mod)

  # output plots (default = HTML with PNG files organized in folders)
  plot_wham_output(mod=mod, dir.main=file.path(getwd(),df.mods$Model[m]), out.type='html')

  # save model
  saveRDS(mod, file=paste0(df.mods$Model[m],".rds"))
}

# collect fit models into a list
mod.list <- paste0(df.mods$Model,".rds")
mods <- lapply(mod.list, readRDS)

# calculate AIC and Mohn's rho
df.aic <- compare_wham_models(mods, sort=FALSE)$tab
df.mods <- cbind(df.mods, df.aic)

# make results table prettier
rownames(df.mods) <- NULL
df.mods$Recruitment <- dplyr::recode(df.mods$Recruitment, `2`='Random', `3`='Bev-Holt', `4`='Ricker')
df.mods$Ecov_how <- dplyr::recode(df.mods$Ecov_how, `0`='---',`1`='Controlling', `2`='Limiting', `4`='Masking')
df.mods <- df.mods[order(df.mods$dAIC),]

# look at results table
df.mods

# save results table
save("df.mods", file="ex2_res.RData")
