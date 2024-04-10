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

# CURRENTLY: selectivity$initial_pars is different between ex4_selectivity and test_ex04_selectivity

is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo
library(ggplot2)
library(tidyr)
library(viridis)
#by default do not perform bias-correction
if(!exists("basic_info")) basic_info <- NULL

# -------------------------------------------------------------------------
# 1. Setup and load data
# -------------------------------------------------------------------------

# create directory for analysis, e.g.
# write.dir <- "/path/to/save/ex2" on linux/mac
if(!exists("write.dir")) write.dir <- tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

# Load data
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))

# --------------------------------------------------------------------------
# 2. Specify selectivity model options
# --------------------------------------------------------------------------
# We are going to run 9 models that differ only in their selectivity options:
# m1-m3 logistic, m4-m8 age-specific
sel_model <- c(rep("logistic",3), rep("age-specific",5))

# time-varying options for each of 3 blocks (b1 = fleet, b2-3 = indices)
sel_re <- list(c("none","none","none"), # m1-m4 logistic
				c("iid","none","none"),
#				c("ar1","none","none"), #can't get a converged model for logistic with two re and two fe for
				c("2dar1","none","none"),
				c("none","none","none"), # m4-m8 age-specific
				c("iid","none","none"),
				c("ar1","none","none"), #will map all age-specific (mean) parameters to one estimated value
				c("ar1_y","none","none"),
				c("2dar1","none","none"))
n.mods <- length(sel_re)

# summary data frame
df.mods <- data.frame(Model=paste0("m",1:n.mods), 
	Selectivity=sel_model, # Selectivity model (same for all blocks)
	Block_1_re=sapply(sel_re, function(x) x[[1]])) # Block 1 random effects
rownames(df.mods) <- NULL

df.mods

#initial pars for logistic or age-specfic selectivity
initial_pars <- c(rep(list(list(c(2,0.3),c(2,0.3),c(2,0.3))),3), rep(list(list(c(0.1,0.5,0.5,1,1,1),c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,1,1,1,1))),5))
#which pars to fix for age-specific selectivity
fix_pars <- c(list(NULL,NULL,NULL), rep(list(list(4:6,4:5,3:6)),5))

# see currently specified selectivity options, from asap3 file:
asap3[[1]]$dat$sel_block_assign # 1 fleet, all years assigned to block 1
# by default each index gets its own selectivity block (here, blocks 2 and 3)

asap3[[1]]$dat$sel_block_option # fleet selectivity (1 block), 2 = logistic
asap3[[1]]$dat$index_sel_option # index selectivity (2 blocks), 2 = logistic
asap3[[1]]$dat$sel_ini # fleet sel initial values (col1), estimation phase (-1 = fix)
asap3[[1]]$dat$index_sel_ini # index sel initial values (col1), estimation phase (-1 = fix)

mods <- list()
# -----------------------------------------------------------------------
# 3. Setup and run models
# -----------------------------------------------------------------------
for(m in 1:n.mods){
	# overwrite initial parameter values in ASAP data file (ex1_SNEMAYT.dat)
	input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2,
		selectivity=list(model=rep(sel_model[m],3), re=sel_re[[m]], initial_pars=initial_pars[[m]], fix_pars = fix_pars[[m]]),
		NAA_re = list(sigma='rec+1',cor='iid'),
		age_comp = "logistic-normal-miss0", basic_info = basic_info) # logistic normal, treat 0 obs as missing

	# fit model
	mods[[m]] <- fit_wham(input, do.check=T, do.osa=F, do.retro=F) 
}

for(m in 1:length(mods)) saveRDS(mod[[m]], file=paste0("m",m,".rds"))

# -----------------------------------------------------------------------
# 4. Model convergence and comparison
# -----------------------------------------------------------------------

# check which models converged
vign4_conv <- lapply(mods, function(x) capture.output(check_convergence(x)))
for(m in 1:n.mods) cat(paste0("Model ",m,":"), vign4_conv[[m]], "", sep='\n')

# plot output for converged models
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
pdHess <- as.logical(ok_sdrep)
conv <- sapply(mods, function(x) x$opt$convergence == 0) # 0 means opt converged
conv_mods <- (1:n.mods)[pdHess] 
for(m in conv_mods){
	plot_wham_output(mod=mods[[m]], dir.main=file.path(write.dir,paste0("m",m)))
}

# compare models using AIC
df.aic <- as.data.frame(compare_wham_models(mods, table.opts=list(sort=FALSE, calc.rho=FALSE))$tab)
df.aic[!pdHess,] = NA
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(data.frame(Model=paste0("m",1:n.mods), Selectivity=sel_model,
	"Block1_re"=sapply(sel_re, function(x) x[[1]]),
	"opt_converged"= ifelse(conv, "Yes", "No"),
	"pd_hessian"= ifelse(pdHess, "Yes", "No"),
	"NLL"=sapply(mods, function(x) round(x$opt$objective,3)),
	"Runtime"=sapply(mods, function(x) x$runtime)), df.aic)
rownames(df.mods) <- NULL
df.mods

# save results table
write.csv(df.mods, file="ex4_table.csv",quote=F, row.names=F)

# plot the models estimates of selectivity-at-age for block 1 (fleet).
# selAA block 1 plots
selAA <- lapply(mods, function(x) x$report()$selAA[[1]])
sel_mod <- factor(c("Age-specific","Logistic")[sapply(mods, function(x) x$env$data$selblock_models[1])], levels=c("Logistic","Age-specific"))
sel_cor <- factor(c("None","IID","AR1","AR1_y","2D AR1")[sapply(mods, function(x) x$env$data$selblock_models_re[1])], levels=c("None","IID","AR1","AR1_y","2D AR1"))
df.selAA <- data.frame(matrix(NA, nrow=0, ncol=11))
colnames(df.selAA) <- c(paste0("Age_",1:6),"Year","Model","conv","sel_mod","sel_cor")
for(m in 1:n.mods){
	df <- as.data.frame(selAA[[m]])
	df$Year <- input$years
	colnames(df) <- c(paste0("Age_",1:6),"Year")
	df$Model <- m
	df$conv <- ifelse(df.mods$pd_hessian[m]=="Yes",1,0)
	df$sel_mod = sel_mod[m]
	df$sel_cor = sel_cor[m]
	df.selAA <- rbind(df.selAA, df)
}
df <- df.selAA %>% pivot_longer(-c(Year,Model,conv,sel_mod,sel_cor),
	names_to = "Age",
	names_prefix = "Age_",
	names_transform = list(Age = as.integer),
	values_to = "Selectivity")
df$Age <- as.integer(df$Age)
df$sel_mod <- factor(df$sel_mod, levels=c("Logistic","Age-specific"))
df$sel_cor <- factor(df$sel_cor, levels=c("None","IID","AR1","AR1_y","2D AR1"))
df$conv = factor(df$conv)

png(file = "selAA.png", width = 8, height = 7.5, res = 200, units='in')
print(ggplot(df, aes(x=Year, y=Age)) +
	geom_tile(aes(fill=Selectivity)) +
	# geom_tile(aes(fill=Selectivity, alpha=conv)) +
	# scale_alpha_discrete(range=c(0.4,1), guide=FALSE) +
	geom_label(aes(x=Year, y=Age, label=lab), size=5, alpha=1, #fontface = "bold",
	  data=data.frame(Year=1975.5, Age=5.8, lab=paste0("m",1:length(mods)), sel_mod=sel_mod, sel_cor=sel_cor)) +	
	facet_grid(rows=vars(sel_cor), cols=vars(sel_mod)) +
	scale_fill_viridis() +
	theme_bw() +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)))
dev.off()
