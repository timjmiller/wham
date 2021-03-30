# WHAM example 8: Read ASAP3 results and compare to fit WHAM models

# source("/home/bstock/Documents/wham/inst/example_scripts/ex8_compare_asap.R")

library(tidyverse)
devtools::load_all("~/Documents/wham/") # test locally (not yet pushed to github)
# devtools::install_github("timjmiller/wham", ref="devel", dependencies=TRUE)
# devtools::install_github("timjmiller/wham", dependencies=TRUE)

# create directory for analysis, e.g.
main.dir <- file.path(getwd(),"sandbox","wham_testing",paste0("runall-",format(Sys.Date(), "%Y%m%d")))
if(!dir.exists(main.dir)) dir.create(main.dir)
write.dir <- file.path(main.dir,"ex8")
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

# Georges Bank haddock BASE_5C.DAT from Liz
asap3 <- read_asap3_dat("/home/bstock/Documents/wham/inst/extdata/BASE_3/BASE_3.DAT")

# Define 4 WHAM models
df.mods <- data.frame(naa_sig=c('none','rec','rec+1','rec+1'), naa_cor=c('none','iid','iid','2dar1'))
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

for(m in 1:n.mods){
  NAA_list <- list(cor=df.mods[m,"naa_cor"], sigma=df.mods[m,"naa_sig"])
  if(NAA_list$sigma == 'none') NAA_list = NULL

  input <- prepare_wham_input(asap3, recruit_model = 2, # Phase for Steepness = -2, Initial Steepness = 1
                              model_name = df.mods$Model[m],                         
                              NAA_re = NAA_list)   

  # Fit model
  mod <- fit_wham(input, do.retro=T, do.osa=F)

  # Save model
  saveRDS(mod, file=file.path(write.dir, paste0(df.mods$Model[m],".rds")))
}

mod.list <- file.path(write.dir,paste0("m",1:n.mods,".rds"))
mods <- lapply(mod.list, readRDS)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- sapply(mods, function(x) x$opt$convergence == 0) # 0 means opt converged
df.mods$pdHess <- as.logical(ok_sdrep)
conv_mods <- (1:n.mods)[df.mods$pdHess] 
for(m in conv_mods){
  plot_wham_output(mod=mods[[m]], out.type='pdf', dir.main=file.path(write.dir,paste0("m",m)))
}

base <- read_asap3_fit(wd="/home/bstock/Documents/wham/inst/extdata/BASE_3", asap.name="BASE_3")

