# WHAM example 4: Time-varying selectivity
# as in example 1
#   no environmental covariate
#   2 indices
#   fit to 1973-2016 data
#   age compositions = 7, logistic normal don't pool zero obs (ex 2 used 5, logistic normal pool` zero obs)
#   selectivity = age-specific
# as in example 2
#   selectivity = logistic

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(tidyverse)
library(viridis)

# -------------------------------------------------------------------------
# 1. Setup and load data
# -------------------------------------------------------------------------

# create directory for analysis, e.g.
# write.dir <- "/path/to/save/ex2" on linux/mac
if(!exists("write.dir")) write.dir = getwd()
dir.create(write.dir)
setwd(write.dir)

# copy `ex1_SNEMAYT.dat` to our analysis directory
wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex1_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)

# Read the ASAP3 .dat file into R and convert to input list for wham:
asap3 <- read_asap3_dat("ex1_SNEMAYT.dat")

# --------------------------------------------------------------------------
# 2. Specify selectivity model options
# --------------------------------------------------------------------------
# We are going to run 10 models that differ only in their selectivity options:
# m1-m5 logistic, m6-m10 age-specific
sel_model <- c(rep("logistic",5), rep("age-specific",5))

# time-varying options for each of 3 blocks (b1 = fleet, b2-3 = indices)
sel_re <- list(c("none","none","none"), # m1-m5 logistic
				c("iid","none","none"),
				c("ar1","none","none"),
				c("ar1_y","none","none"),
				c("2dar1","none","none"),
				c("none","none","none"), # m6-m10 age-specific
				c("iid","none","none"),
				c("ar1","none","none"),
				c("ar1_y","none","none"),
				c("2dar1","none","none"))
n.mods <- length(sel_re)

# summary data frame
df.mods <- data.frame(Model=paste0("m",1:n.mods), Selectivity=sel_model, 
							"Block_1"=sapply(sel_re, function(x) x[[1]]),
							"Block_2"=sapply(sel_re, function(x) x[[2]]),
							"Block_3"=sapply(sel_re, function(x) x[[3]]))
rownames(df.mods) <- NULL

# look at model table
df.mods

# see currently specified selectivity options, from asap3 file:
asap3$dat$sel_block_assign # 1 fleet, all years assigned to block 1
# by default each index gets its own selectivity block (here, blocks 2 and 3)

asap3$dat$sel_block_option # fleet selectivity (1 block), 2 = logistic
asap3$dat$index_sel_option # index selectivity (2 blocks), 2 = logistic
asap3$dat$sel_ini # fleet sel initial values (col1), estimation phase (-1 = fix)
asap3$dat$index_sel_ini # index sel initial values (col1), estimation phase (-1 = fix)

# -----------------------------------------------------------------------
# 3. Setup and run models
# -----------------------------------------------------------------------
inv.logit <- function(x) exp(x)/(1+exp(x))
mods <- vector("list",n.mods) # store models in a list
selAA <- vector("list",n.mods) # save selectivity-at-age for block 1 for each model
for(m in 1:n.mods){
	if(sel_model[m] == "logistic"){ # logistic selectivity
		# overwrite initial parameter values in ASAP data file (ex1_SNEMAYT.dat)
		input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2, 
					selectivity=list(model=rep("logistic",3), re=sel_re[[m]], initial_pars=list(c(inv.logit(-0.67935549),0.2),c(2,0.2),c(2,0.2))))
		input$par$sel_repars[1,1] <- -1.3
	} else { # age-specific selectivity
		# fix ages 1,4,5 / 4 / 2
		input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2, 
					selectivity=list(model=rep("age-specific",3), re=sel_re[[m]], initial_pars=list(c(inv.logit(-4),0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,0.5,0.5,0.5,0.5)), fix_pars=list(c(1,4,5),4,2)))
		input$par$sel_repars[1,1] <- -0.4
	}

	# overwrite age comp model (all models use logistic normal)
	input$data$age_comp_model_indices = rep(7, input$data$n_indices)
	input$data$age_comp_model_fleets = rep(7, input$data$n_fleets)
	input$data$n_age_comp_pars_indices = rep(1, input$data$n_indices)
	input$data$n_age_comp_pars_fleets = rep(1, input$data$n_fleets)
	input$par$index_paa_pars = rep(0, input$data$n_indices)
	input$par$catch_paa_pars = rep(0, input$data$n_fleets)
	input$map = input$map[!(names(input$map) %in% c("index_paa_pars", "catch_paa_pars"))]

	# overwrite NAA model (all models use full state space)
	input$data$use_NAA_re = 1
	input$data$random_recruitment = 0
	input$map = input$map[!(names(input$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
	input$map$log_R = factor(rep(NA, length(input$par$log_R)))
	input$random = c(input$random, "log_NAA")
	
	# fit model
	mods[[m]] <- fit_wham(input, do.check=T, do.osa=F, do.proj=F, do.retro=F) 
	saveRDS(mods[[m]], file=paste0("m",m,".rds"))

	# save selectivity-at-age for block 1 (fleet)
	selAA[[m]] <- mods[[m]]$report()$selAA[[1]]
}

# -----------------------------------------------------------------------
# 4. Model convergence and comparison
# -----------------------------------------------------------------------
# check models converged
sapply(mods, function(x) check_convergence(x))

# compare models using AIC
df.aic <- compare_wham_models(mods, sort=FALSE, calc.rho=FALSE)$tab
df.mods <- cbind(data.frame(Model=paste0("m",1:n.mods), Selectivity=sel_model, 
							"Block_1"=sapply(sel_re, function(x) x[[1]]),
							"Block_2"=sapply(sel_re, function(x) x[[2]]),
							"Block_3"=sapply(sel_re, function(x) x[[3]]),
							"NLL"=sapply(mods, function(x) round(x$opt$objective,3))), df.aic)
rownames(df.mods) <- NULL
df.mods

# plot the models estimates of selectivity-at-age for block 1 (fleet).
df.selAA <- data.frame(matrix(NA, nrow=0, ncol=8))
colnames(df.selAA) <- c(paste0("Age_",1:6),"Year","Model")
for(m in 1:n.mods){
	df <- as.data.frame(selAA[[m]])
	df$Year <- input$years
	colnames(df) <- c(paste0("Age_",1:6),"Year")
	df$Model <- m
	df.selAA <- rbind(df.selAA, df)
}
df <- df.selAA %>% pivot_longer(-c(Year,Model),
				names_to = "Age", 
				names_prefix = "Age_",
				names_ptypes = list(Age = integer()),
				values_to = "Selectivity")
df$sel_model <- factor(rep(c("Logistic","Age-specific"), each=dim(df)[1]/2), levels=c("Logistic","Age-specific"))
df$sel_re <- factor(c(rep(c("None","IID","AR1","AR1_y","2D AR1"), each=dim(df)[1]/n.mods), rep(c("None","IID","AR1","AR1_y","2D AR1"), each=dim(df)[1]/n.mods)), levels=c("None","IID","AR1","AR1_y","2D AR1"))

print(ggplot(df, aes(x=Year, y=Age, fill=Selectivity)) + 
	geom_tile() +
	theme_bw() + 
	facet_grid(rows=vars(sel_re), cols=vars(sel_model)) +
	scale_fill_viridis())
