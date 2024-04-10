# WHAM example 8: Read ASAP3 results and compare to fit WHAM models

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo
#by default do not perform bias-correction
if(!exists("basic_info")) basic_info <- NULL

library(tidyverse)

# create directory for analysis, e.g.
if(!exists("write.dir")) write.dir <- tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

# Georges Bank haddock BASE_5C.DAT
path_to_examples <- system.file("extdata", package="wham")

# Georges Bank haddock BASE_5C.DAT
asap3 <- read_asap3_dat(file.path(path_to_examples,"BASE_3", "BASE_3.DAT"))

# Define WHAM models
#   Model naa_sig naa_cor
#      m1    none    none  fixed effect recruitment
#      m2     rec     iid  random effect recruitment
#      m3   rec+1     iid  all numbers-at-age are random effects, full state-space model
df.mods <- data.frame(naa_sig=c('none','rec','rec+1'), naa_cor=c('none','iid','iid'))
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

mods <- list()
for(m in 1:n.mods){
  NAA_list <- list(cor=df.mods[m,"naa_cor"], sigma=df.mods[m,"naa_sig"])
  if(NAA_list$sigma == 'none') NAA_list = NULL

  input <- prepare_wham_input(asap3, recruit_model = 2, # match asap model, which does not estimate stock-recruitment relationship
                              model_name = df.mods$Model[m],                         
                              NAA_re = NAA_list, basic_info = basic_info)   

  mods[[m]] <- fit_wham(input, do.osa=F)
  saveRDS(mods[[m]], file=file.path(write.dir, paste0(df.mods$Model[m],".rds")))
}

ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- sapply(mods, function(x) x$opt$convergence == 0) # 0 means opt converged
df.mods$pdHess <- as.logical(ok_sdrep)
conv_mods <- (1:n.mods)[df.mods$pdHess] 
for(m in conv_mods){
  plot_wham_output(mod=mods[[m]], out.type='pdf', dir.main=file.path(write.dir,paste0("m",m)))
}

# get output from ASAP model run
base <- read_asap3_fit(wd=file.path(path_to_examples,"BASE_3"), asap.name="BASE_3")
mods <- c(list(base),mods)
names(mods) <- c("ASAP",paste0("WHAM-",df.mods$Model))

# compare models
res <- compare_wham_models(mods, fdir=write.dir, plot.opts=list(kobe.prob=FALSE))
saveRDS(res, file=file.path(write.dir,"res.rds"))

# aic table
round(res$tab,2)

# --------------------------------------------------------------------------
# lots of options, see ?compare_wham_models
# only get table, not plots
res <- compare_wham_models(mods, fdir=write.dir, do.plot=F)

# only get plots, not table
res <- compare_wham_models(mods, fdir=write.dir, do.table=F)

# plot.opts list of options for plots:
#   $out.type - character, either \code{'pdf'} or \code{'png'} (default = \code{'png'} because I am not sure \code{system('pdftk')} will work across platforms.)}
#   $ci - vector of T/F, length = 1 (applied to all models) or number of models}
#   $years - vector, which years to plot? Default = all (model and projection years).}
#   $which - vector, which plots to make? Default = all. See details.}
#   $relative.to - character, name of "base" model to plot differences relative to.}
#   $alpha - scalar, (1-alpha)\% confidence intervals will be plotted. Default = 0.05 for 95\% CI.}
#   $ages.lab - vector, overwrite model age labels.}
#   $kobe.yr - integer, which year to use in Kobe plot (relative status). Default = terminal model year.}
#   $M.age - integer, which age to use in M time-series plot. Default = `data$which_F_age` (age of F to use for full total F).}
#   $return.ggplot - T/F, return a list of ggplot2 objects for later modification? Default = TRUE.}
#   $kobe.prob - T/F, print probabilities for each model in each quadrant of Kobe plot? Default = TRUE.}

# 10 plots are produced
#   1) 3-panel of SSB (spawning stock biomass), F (fully-selected fishing mortality), and Recruitment
#   2) CV (coefficient of variation) for SSB, F, and Recruitment
#   3) Fleet selectivity (by block, averaged across years)
#   4) Index selectivity (by block, averaged across years)
#   5) Selectivity tile (fleets + indices, useful for time-varying random effects)
#   6) M time series (natural mortality, can specify which age with plot.opts$M.age)
#   7) M tile (useful for time-varying random effects)
#   8) 3-panel of F X\% SPR, SSB at F_X\%SPR, and yield at F_X\%SPR
#   9) 2-panel of relative status (SSB / SSB at F_X\%SPR and F / F_X\%SPR)
#   10) Kobe status (relative SSB vs. relative F)

# plots are saved as png by default, can be pdf
res <- compare_wham_models(mods, fdir=write.dir, plot.opts=list(out.type='pdf'))

# which = 9 (only plot relative status)
# years = 1980-2018
compare_wham_models(mods, do.table=F, plot.opts=list(years=1980:2018, which=9))
ggsave(file.path(write.dir,"which9_zoom.png"), device='png', width=6.5, height=5.5, units='in')

# which = 1 (SSB, F, recruitment)
# ci = FALSE (remove confidence intervals for all models, can also choose a subset)
compare_wham_models(mods, fdir=write.dir, do.table=F, plot.opts=list(ci=FALSE, which=1))
ggsave(file.path(write.dir,"which1_noCI.png"), device='png', width=6.5, height=5.5, units='in')

# which = 2 (CV of SSB, F, recruitment)
# relative to ASAP
compare_wham_models(mods, fdir=write.dir, do.table=F, plot.opts=list(ci=FALSE, relative.to="ASAP", which=2))
ggsave(file.path(write.dir,"which2_relative.png"), device='png', width=6.5, height=5.5, units='in')

# which = 10 (kobe plot)
# kobe.yr = 2010 (instead of terminal year, 2018)
# kobe.prob = F (don't print probabilities)
compare_wham_models(mods, fdir=write.dir, do.table=F, plot.opts=list(which=10, kobe.yr=2010, kobe.prob=F))

# res$g holds the ggplot objects so you can modify later
# relative status plot with different fill and color scales
res$g[[9]] + scale_colour_brewer(palette="Set1") + scale_fill_brewer(palette="Set1")
ggsave(file.path(write.dir,"which9_colorchange.png"), device='png', width=6.5, height=5.5, units='in')

# F, SSB, Recruitment of WHAM models relative to ASAP
# note: even if only making one plot, res$g is still length(10) and we want res$g[[1]]
# zoom in on 1980-2018, remove CI
res <- compare_wham_models(mods, do.table=F, plot.opts=list(years=1980:2018, ci=FALSE, relative.to="ASAP", which=1))
cols <- c("black", RColorBrewer::brewer.pal(n = 3, name = "Set1"))
res$g[[1]] + scale_colour_manual(values=cols)
ggsave(file.path(write.dir,"which1_relative_colorchange.png"), device='png', width=5, height=5.5, units='in')

# any aesthetics that weren't in the original plot are more complicated
# for example if we want to make the base model line dashed, linetype was not in original plot
res$g[[1]]$mapping$linetype = quote(Model)
res$g[[1]]$labels$linetype = "Model"
ltys <- c(2,1,1,1)
res$g[[1]] + scale_colour_manual(values=cols) + scale_linetype_manual(values=ltys)
ggsave(file.path(write.dir,"which1_relative_colorchange_linetype.png"), device='png', width=5, height=5.5, units='in')

# change labels on index selectivity plot facets
res <- compare_wham_models(mods, do.table=F, plot.opts=list(which=4))
index_names <- as_labeller(c(`Block 4` = "NEFSC - Spring", `Block 5` = "NEFSC - Fall",`Block 6` = "DFO", `Block 7` = "NEFSC - Spring41"))
res$g[[4]] + facet_wrap(vars(Block), ncol=1, strip.position = 'right', labeller = index_names)
ggsave(file.path(write.dir,"which4_labels.png"), device='png', width=4.5, height=5.5, units='in')
