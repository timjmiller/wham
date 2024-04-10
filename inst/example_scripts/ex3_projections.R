# WHAM example 3: projections
# Replicate Miller et al 2016 results
#   adds environmental covariate (CPI, treated as rw)
#   uses 5 indices (ex 1: only 2 indices)
#   only fit to 1973-2011 data (ex 1: 1973-2016)
#   age compositions = 5, logistic normal pool zero obs (ex 1: 7, logistic normal missing zero obs)
#   selectivity = logistic (ex 1: age-specific)

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo
wham.dir <- find.package("wham")
#by default do not perform bias-correction
if(!exists("basic_info")) basic_info <- NULL


# create directory for analysis, e.g.
if(!exists("write.dir")) write.dir <- tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

# ---------------------------------------------------------
# Load data
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples, "CPI.csv"), header=T)

# specify model: AR1 CPI, limiting effect on Bev-Holt
env <- list(
  label = "CPI",
  mean = as.matrix(env.dat$CPI), # CPI observations
  logsigma = as.matrix(log(env.dat$CPI_sigma)), # CPI standard error is given/fixed as data
  year = env.dat$Year,
  use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
  process_model = "ar1", # fit CPI as AR1 process
  recruitment_how = matrix("limiting-lag-1-linear")) # limiting (carrying capacity), CPI in year t affects recruitment in year t+1

input <- prepare_wham_input(asap3, recruit_model = 3,
                            model_name = "Ex 3: Projections",
                            ecov = env,
                            NAA_re = list(sigma="rec+1", cor="iid"),
                            age_comp = "logistic-normal-pool0", basic_info = basic_info) # logistic normal pool 0 obs

# selectivity = logistic, not age-specific
#   2 pars per block instead of n.ages
#   sel pars of indices 4/5 fixed at 1.5, 0.1 (neg phase in .dat file)
input$par$logit_selpars[1:4,7:8] <- 0 # original code started selpars at 0 (last 2 rows are fixed)
#input$data$FMSY_init[] = 0.2
# ---------------------------------------------------------
## Fit model without projections

mod <- fit_wham(input) # default do.proj=FALSE
saveRDS(mod, file="m5.rds")
#mod <- readRDS("m5.rds")

# Add projections to previously fit model
mod_proj <- list()

# default settings: 3 years, use last F, continue ecov
mod_proj[[1]] <- project_wham(mod, proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE)

# 5 years, use last F, average ecov 1992-1996
mod_proj[[2]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=1992:1996, proj.ecov=NULL), save.sdrep=FALSE)

# 5 years, use last F, use last ecov
mod_proj[[3]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=TRUE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE)

# 5 years, use last F, specify high CPI ~ 0.5
mod_proj[[4]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=matrix(c(0.5,0.7,0.4,0.5,0.55),ncol=1)), save.sdrep=FALSE)

# 5 years, use last F, specify low CPI ~ -1.5
mod_proj[[5]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=matrix(c(-1.6,-1.3,-1,-1.2,-1.25),ncol=1)), save.sdrep=FALSE)

# specify catch, 5 years
mod_proj[[6]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=c(10, 2000, 1000, 3000, 20), avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE)

# specify F, 5 years
mod_proj[[7]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=c(0.001, 1, 0.5, .1, .2), proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE)

# use FXSPR (avg.yrs defaults to last 5 years, 2007-2011), 5 years
mod_proj[[8]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE)

# use avg F (avg.yrs defaults to last 5 years, 2007-2011), 3 years
mod_proj[[9]] <- project_wham(mod, proj.opts=list(n.yrs=3, use.last.F=FALSE, use.avg.F=TRUE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE)

# use avg F 1992-1996, 10 years
mod_proj[[10]] <- project_wham(mod, proj.opts=list(n.yrs=10, use.last.F=FALSE, use.avg.F=TRUE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=1992:1996,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE)

# use FMSY (avg.yrs defaults to last 5 years, 2007-2011), 5 years
mod_proj[[11]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FMSY=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE)


saveRDS(mod_proj, file="m5_proj.rds")
# mod_proj <- readRDS("m5_proj.rds")

#  check marginal nll is the same
nll_proj <-  sapply(mod_proj, function(x) x$opt$obj)
mod$opt$obj
round(nll_proj - mod$opt$obj, 6) # difference between original and projected models' NLL

# plot results
for(m in 1:length(mod_proj)){
  plot_wham_output(mod_proj[[m]], dir.main=file.path(write.dir,paste0("proj_",m)))
}

# to more easily compare plots, copy to folders organized by plot instead of by model
# plots <- c("Ecov_1","F_byfleet","SSB_at_age","SSB_F_trend","SSB_Rec_time","Kobe_status")
# dirs <- file.path(write.dir,plots)
# lapply(as.list(dirs), FUN=dir.create)
# for(m in 1:length(mod_proj)){
#   for(i in 1:(length(plots)-1)){
#      file.copy(from=file.path(write.dir,paste0("proj_",m),"plots_png","results",paste0(plots[i],".png")),
#                to=file.path(dirs[i],paste0(plots[i],"_proj_",m,".png")))
#   }
#   file.copy(from=file.path(write.dir,paste0("proj_",m),"plots_png","ref_points","Kobe_status.png"),
#             to=file.path(dirs[i+1],paste0(plots[i+1],"_proj_",m,".png")))
# }
