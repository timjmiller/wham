# WHAM example 9: Retrospective predictions

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
# devtools::load_all("~/Documents/wham/")
is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo
#by default do not perform bias-correction
if(!exists("basic_info")) basic_info <- NULL

library(tidyverse)
library(viridis)
library(ggplot2)

# create directory for analysis, e.g.
if(!exists("write.dir")) write.dir <- tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

# Load data
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"CPI.csv"), header=T)

# 2 models, with and without CPI effect on recruitment (both fit CPI data to compare AIC)
# Model  Recruit_mod  Ecov_mod     Ecov_how
#    m1     Bev-Holt       ar1        ---
#    m2     Bev-Holt       ar1     Limiting
df.mods <- data.frame(Model = c("m1","m2"), Ecov_how = c(0,2), 
	R_how = c("none", "limiting-lag-1-linear"), stringsAsFactors=FALSE)

n.mods <- dim(df.mods)[1]
df.mods

# specify CPI model
env <- list(
  label = "CPI",
  mean = as.matrix(env.dat$CPI), # CPI observations
  logsigma = as.matrix(log(env.dat$CPI_sigma)), # CPI standard error is given/fixed as data
  year = env.dat$Year,
  use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
  process_model = 'ar1')#, # "rw" or "ar1"

mods <- list()
for(m in 1:n.mods){
	env$recruitment_how = matrix(df.mods$R_how[m],1,1)
	input <- prepare_wham_input(asap3, recruit_model = 3,
	                            model_name = df.mods$Model[m],
	                            ecov = env,
	                            NAA_re = list(sigma="rec+1", cor="iid"),
	                            age_comp = "logistic-normal-pool0", basic_info = basic_info) # logistic normal pool 0 obs
	input$par$logit_selpars[1:4,7:8] <- 0

	# fit model
	mods[[m]] <- fit_wham(input, do.retro=F, do.osa=F)
	saveRDS(mods[[m]], file=paste0(df.mods$Model[m],".rds"))
	plot_wham_output(mod=mods[[m]], dir.main=file.path(write.dir,df.mods$Model[m]), out.type='png')
}

ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- sapply(mods, function(x) x$opt$convergence == 0) # 0 means opt converged
df.mods$pdHess <- as.logical(ok_sdrep)
df.mods

names(mods) <- c("m1 (no CPI)","m2 (CPI)")
res <- compare_wham_models(mods, fdir=write.dir, do.table=F, plot.opts=list(kobe.prob=FALSE))

n.yrs.peel <- 15
n.yrs.proj <- 3
for(m in 1:n.mods){
	mods[[m]]$peels <- retro(mods[[m]], ran = unique(names(mods[[m]]$env$par[mods[[m]]$env$random])), n.peels=n.yrs.peel, save.input = T)
	for(p in 3:n.yrs.peel) {
		mods[[m]]$peels[[p]] <- project_wham(mods[[m]]$peels[[p]], proj.opts = list(n.yrs = n.yrs.proj, proj.F=rep(0.001,n.yrs.proj)), check.version = FALSE)
	}
}
# plot retrospective predictions of recruitment
plot_retro_pred_R <- function(mods, peels=3:n.yrs.peel, n.yrs.proj=n.yrs.proj){
	df <- data.frame(matrix(NA, nrow=0, ncol=6))
	colnames(df) <- c("Year","Model","Peel","Recruitment","termyr")
	for(m in 1:length(mods)){
		for(p in peels){
			tmp <- suppressWarnings(read_wham_fit(mods[[m]]$peels[[p]]))
			df <- rbind(df, data.frame(Year=tail(tmp$years_full, n.yrs.proj+1),
										Model=names(mods)[m],
										Peel=p,
										Recruitment=tail(exp(tmp$log_NAA_rep$est[1,1,,1]), n.yrs.proj+1),
										termyr=c(1,rep(0,n.yrs.proj))))
		}
		# get full model fit, "peel 0"
		tmp <- suppressWarnings(read_wham_fit(mods[[m]]))
		df <- rbind(df, data.frame(Year=tmp$years,
									Model=names(mods)[m],
									Peel=0,
									Recruitment=exp(tmp$log_NAA_rep$est[1,1,1:length(tmp$years),1]),
									termyr=0))
	}

	df <- filter(df, Year > 1990)
	df$Model <- factor(df$Model, levels=names(mods), labels=names(mods))
	df$Year <- as.integer(df$Year) 
	df$Peel <- factor(df$Peel)
	dfpts <- filter(df, termyr==1)

	cols <- c("black", viridis_pal(option="plasma")(length(levels(df$Peel))-1))
	g <- ggplot(df, aes(x=Year, y=Recruitment, linetype=Model, color=Peel, fill=Peel, group=interaction(Model,Peel))) + 
	      geom_line(linewidth=1) +
	      geom_point(data=dfpts, color='black', size=2, pch=21) +
	      scale_x_continuous(expand=c(0.01,0.01)) + # breaks=scales::breaks_extended(5)
	      scale_colour_manual(values=cols) +
	      scale_fill_manual(values=cols) +
	      guides(color = "none", fill="none") +
	      scale_y_continuous(expand=c(0.01,0.01), limits = c(0,NA), labels=fancy_scientific) +
	      theme_bw() +
	      theme(legend.position=c(.9,.9), legend.box.margin = margin(0,0,0,0), legend.margin = margin(0,0,0,0))
	return(g)
}

plot_retro_pred_R(mods, peels=3:n.yrs.peel, n.yrs.proj=n.yrs.proj)
ggsave(file.path(write.dir,"retro_pred_R.png"), device='png', width=8, height=5, units='in')

