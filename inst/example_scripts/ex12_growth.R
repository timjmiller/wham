## ---- include = FALSE---------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
wham.dir <- find.package("wham")
knitr::opts_knit$set(root.dir = file.path(wham.dir,"extdata"))


## ----message=FALSE------------------------------------------------------------------------------------------------------
library(wham)


## ----eval=FALSE---------------------------------------------------------------------------------------------------------
## wham.dir <- find.package("wham")
## file.path(wham.dir, "example_scripts")


## ---- eval=FALSE--------------------------------------------------------------------------------------------------------
## write.dir <- "choose/where/to/save/output" # otherwise will be saved in working directory
## source(file.path(wham.dir, "example_scripts", "ex12_growth.R"))


## ---- eval=FALSE--------------------------------------------------------------------------------------------------------
## # choose a location to save output, otherwise will be saved in working directory
## write.dir <- "choose/where/to/save/output" # need to change
## dir.create(write.dir)
## setwd(write.dir)


## ----eval=FALSE---------------------------------------------------------------------------------------------------------
## wham.dir <- find.package("wham")
## file.copy(from=file.path(wham.dir,"extdata","ex12_data_file.rds"), to=write.dir, overwrite=FALSE)
## file.copy(from=file.path(wham.dir,"extdata","ex12_SS_report.rds"), to=write.dir, overwrite=FALSE)
## file.copy(from=file.path(wham.dir,"extdata","PDO_annual_1970_2019.csv"), to=write.dir, overwrite=FALSE)


## -----------------------------------------------------------------------------------------------------------------------
list.files()


## -----------------------------------------------------------------------------------------------------------------------
data_file <- readRDS("ex12_data_file.rds")
SS_report <- readRDS("ex12_SS_report.rds")
envindex = read.csv("PDO_annual_1970_2019.csv")


## -----------------------------------------------------------------------------------------------------------------------
min_year = SS_report$startyr # first year
max_year = SS_report$endyr # last year
n_ages = length(1:max(SS_report$agebins)) # number of ages
fish_len = SS_report$lbins # length bins
n_years = length(min_year:max_year) # number of years
NAA_SS = SS_report$natage[SS_report$natage$`Beg/Mid` == 'B' & SS_report$natage$Yr >= min_year & SS_report$natage$Yr <= max_year, 14:ncol(SS_report$natage)] # abundance at age from SS
LAA_SS = SS_report$growthseries[SS_report$growthseries$Yr >= min_year & SS_report$growthseries$Yr <= max_year, 6:ncol(SS_report$growthseries)] # mean length at age from SS
LWpars = c(SS_report$Growth_Parameters$WtLen1, SS_report$Growth_Parameters$WtLen2) # LW parameters
GWpars = c(SS_report$Growth_Parameters$K, SS_report$Growth_Parameters$Linf, 
           SS_report$endgrowth$Len_Beg[SS_report$endgrowth$Real_Age == 1], # This is L1
           SS_report$endgrowth$SD_Beg[2], SS_report$endgrowth$SD_Beg[11]) # growth parameters (vB function)
Q_pars = c(exp(SS_report$parameters$Value[SS_report$parameters$Label == "LnQ_base_Fleet2(2)"])) # Q parameter
# Selectivity parameters:
int_selPos = grep(pattern =  "Size_inflection", x = SS_report$parameters$Label)[1]
SelecParams = SS_report$parameters[int_selPos:nrow(SS_report$parameters), ]
selpars1 = SelecParams$Value[1:2] # logistic
selpars2 = SelecParams$Value[3:4] # logistic


## -----------------------------------------------------------------------------------------------------------------------
wham_data = list()
wham_data$ages = 1:n_ages
wham_data$lengths = fish_len
wham_data$years = min_year:max_year

wham_data$n_fleets = 1L
wham_data$agg_catch = as.matrix(data_file$catch[2:(n_years+1), 4])
wham_data$catch_cv = as.matrix(data_file$catch[2:(n_years+1), 5])
# Length information fleet
catch_pal_num = data_file$lencomp[data_file$lencomp$FltSvy == 1, 7:(7+length(wham_data$lengths)-1)]
catch_pal_prop = t(apply(as.matrix(catch_pal_num),1, function(x) x/sum(x)))
wham_data$catch_pal = catch_pal_prop
wham_data$catch_NeffL = as.matrix(as.double(data_file$lencomp[data_file$lencomp$FltSvy == 1, 6]))
wham_data$use_catch_pal = matrix(1, nrow = length(wham_data$years), ncol = wham_data$n_fleets)
wham_data$selblock_pointer_fleets = matrix(1L, ncol = 1, nrow = (max_year - min_year + 1))
wham_data$F = matrix(apply(X = SS_report$fatage[SS_report$fatage$Yr >= min_year &
                                                  SS_report$fatage$Yr <= max_year,
                                                9:ncol(SS_report$fatage)], MARGIN = 1, FUN = max),
                     ncol = 1)

wham_data$n_indices = 1L
wham_data$agg_indices = as.matrix(data_file$CPUE[,4])
# Length information survey
index_pal_num = data_file$lencomp[data_file$lencomp$FltSvy == 2, 7:(7+length(wham_data$lengths)-1)]
index_pal_prop = t(apply(as.matrix(index_pal_num),1, function(x) x/sum(x)))
wham_data$index_pal = array(data = index_pal_prop, dim = c(1, dim(index_pal_prop)[1], dim(index_pal_prop)[2]))
wham_data$index_NeffL =  as.matrix(as.double(data_file$lencomp[data_file$lencomp$FltSvy == 2, 6]))
wham_data$use_index_pal = matrix(1, ncol = 1, nrow = (max_year - min_year + 1))
wham_data$index_cv = as.matrix(data_file$CPUE[,5])
wham_data$units_indices = rep(0L, times = 1) # numbers
wham_data$selblock_pointer_indices = matrix(2L, ncol = 1, nrow = (max_year - min_year + 1))
wham_data$fracyr_indices = matrix(0, ncol = 1, nrow = (max_year - min_year + 1))
wham_data$maturity = matrix(rep(SS_report$endgrowth[2:(n_ages+1),19], times = max_year - min_year + 1),
                            ncol = n_ages, nrow = max_year - min_year + 1, byrow = TRUE) 
wham_data$fracyr_SSB = rep(0, times = max_year - min_year + 1)
# Age comps index:
wham_data$index_paa = array(NA, dim = c(wham_data$n_indices, n_years, n_ages))
this_matrix = data_file$agecomp[,11:(10+n_ages)]
wham_data$index_paa[1,,] = as.matrix(this_matrix/rowSums(this_matrix))
wham_data$index_Neff = matrix(data_file$agecomp$Nsamp,
                              nrow = n_years, ncol = 1)
wham_data$use_index_paa = matrix(1L, nrow = n_years, ncol = 1)
# WAA information 
wham_data$waa_pointer_indices = 1
wham_data$waa_pointer_fleets = 2
wham_data$waa_pointer_totcatch = 2
wham_data$waa_pointer_ssb = 1
wham_data$waa_pointer_jan1 = 1
wham_data$Fbar_ages = 1L:20L
wham_data$percentSPR = 40
wham_data$percentFXSPR = 100
wham_data$percentFMSY = 100
wham_data$XSPR_R_avg_yrs = 1:n_years
wham_data$XSPR_R_opt = 2
wham_data$simulate_period = c(1,0)
wham_data$bias_correct_process = 1
wham_data$bias_correct_observation = 1


## -----------------------------------------------------------------------------------------------------------------------
my_input1 = prepare_wham_input(model_name="whamGrowth1",
                               selectivity=list(model=rep("len-logistic",2), 
                                                re=rep("none",2), 
                                                initial_pars=list(selpars1,selpars2),
                                                fix_pars = list(2, 2),
                                                n_selblocks = 2),
                              M = list(model = 'constant', re = 'none',
                                       initial_means = SS_report$Natural_Mortality[1,5],
                                       est_ages = 1),
                              NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                            recruit_model = 2,
                                            N1_pars = c(NAA_SS[1,1], 0),
                                            recruit_pars = mean(NAA_SS[,1])),
                              growth = list(model = 'vB_classic',
                                            re = c('iid_y', 'none', 'none'),
                                            init_vals = GWpars[1:3],
                                            est_pars = 1:3,
                                            SD_vals = GWpars[4:5], SD_est = 1:2),
                              LW = list(init_vals = LWpars,
                                        re = c('none', 'none')),
                              catchability = list(re = 'none', initial_q = Q_pars, q_lower = 0,
                                                  q_upper = 10, prior_sd = NA),
                              basic_info = wham_data)


## -----------------------------------------------------------------------------------------------------------------------
my_input1$par$log_NAA = as.matrix(log(NAA_SS[-1,])) # set initial rec devs as OM
my_input1$par$log_NAA_sigma = log(SS_report$sigma_R_in) # sigma as in SS
my_input1$map$log_NAA_sigma = factor(NA) # fix sigma
my_input1$map$log_N1_pars = factor(c(1,NA))
my_input1$random = "growth_re"


## -----------------------------------------------------------------------------------------------------------------------
my_model_1 = fit_wham(my_input1, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE)
check_convergence(my_model_1)


## -----------------------------------------------------------------------------------------------------------------------
my_input2 = prepare_wham_input(model_name="whamGrowth2",
                               selectivity=list(model=rep("len-logistic",2), 
                                                re=rep("none",2), 
                                                initial_pars=list(c(40, 5),c(20,3)),
                                                fix_pars = list(NULL, NULL),
                                                n_selblocks = 2),
                               M = list(model = 'constant', re = 'none',
                                        initial_means = SS_report$Natural_Mortality[1,5],
                                        est_ages = 1),
                               NAA_re = list(sigma="rec", cor = 'iid', N1_model = 0,
                                             recruit_model = 2,
                                             N1_pars = as.vector(as.matrix(NAA_SS[1,])),
                                             recruit_pars = mean(NAA_SS[,1])),
                               LW = list(init_vals = LWpars,
                                         re = c('none', 'none')),
                               LAA = list(LAA_vals = as.vector(colMeans(LAA_SS)),
                                          re = c('iid'),
                                          SD_vals = GWpars[4:5], SD_est = 1:2),
                               catchability = list(re = 'none', initial_q = 1, q_lower = 0,
                                                   q_upper = 10, prior_sd = NA),
                               basic_info = wham_data)


## -----------------------------------------------------------------------------------------------------------------------
my_input2$par$log_NAA = as.matrix(log(NAA_SS[-1,])) # set initial rec devs as OM
my_input2$par$log_NAA_sigma = log(SS_report$sigma_R_in) # sigma as in SS
my_input2$map$log_NAA_sigma = factor(NA) # fix sigma
my_input2$map$log_N1_pars = factor(c(1,NA))
my_input2$random = "LAA_re"


## -----------------------------------------------------------------------------------------------------------------------
my_model_2 = fit_wham(my_input2, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE)
check_convergence(my_input2)


## -----------------------------------------------------------------------------------------------------------------------
envindex = envindex[envindex$years >= min_year & envindex$years <= max_year, ]
ecov <- list(
  label = "PDO",
  mean = as.matrix(envindex$index),
  logsigma = 'est_1',
  year = envindex$years,
  use_obs = matrix(1, ncol=1, nrow=dim(envindex)[1]),
  lag = 0,
  ages = list(1:n_ages),
  process_model = 'ar1', # "rw" or "ar1"
  where = 'growth', 
  where_subindex = 1, # K
  how = 1)


## -----------------------------------------------------------------------------------------------------------------------
my_input3 = prepare_wham_input(model_name="whamGrowth1",
                               selectivity=list(model=rep("len-logistic",2), 
                                                re=rep("none",2), 
                                                initial_pars=list(selpars1,selpars2),
                                                fix_pars = list(2, 2),
                                                n_selblocks = 2),
                              M = list(model = 'constant', re = 'none',
                                       initial_means = SS_report$Natural_Mortality[1,5],
                                       est_ages = 1),
                              NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                            recruit_model = 2,
                                            N1_pars = c(NAA_SS[1,1], 0),
                                            recruit_pars = mean(NAA_SS[,1])),
                              growth = list(model = 'vB_classic',
                                            re = c('none', 'none', 'none'),
                                            init_vals = GWpars[1:3],
                                            est_pars = 1:3,
                                            SD_vals = GWpars[4:5], SD_est = 1:2),
                              ecov = ecov,
                              LW = list(init_vals = LWpars,
                                        re = c('none', 'none')),
                              catchability = list(re = 'none', initial_q = Q_pars, q_lower = 0,
                                                  q_upper = 10, prior_sd = NA),
                              basic_info = wham_data)


## -----------------------------------------------------------------------------------------------------------------------
my_input3$par$log_NAA = as.matrix(log(NAA_SS[-1,])) # set initial rec devs as OM
my_input3$par$log_NAA_sigma = log(SS_report$sigma_R_in) # sigma as in SS
my_input3$map$log_NAA_sigma = factor(NA) # fix sigma
my_input3$map$log_N1_pars = factor(c(1,NA))
my_input3$random = "Ecov_re"


## -----------------------------------------------------------------------------------------------------------------------
my_model_3 = fit_wham(my_input3, do.osa = FALSE, do.fit = TRUE, do.retro = FALSE)
check_convergence(my_model_3)


## -----------------------------------------------------------------------------------------------------------------------
LAA_SS = as.matrix(LAA_SS) # true
rownames(LAA_SS) = min_year:(max_year-1)
data_0 = setNames(reshape2::melt(LAA_SS), c('year', 'age', 'LAA'))
data_0$age = data_0$age
data_0$type = 'OM'
# Model 1:
LAA_data_1 = my_model_1$rep$LAA
rownames(LAA_data_1) = min_year:max_year
colnames(LAA_data_1) = 1:n_ages
data_1 = setNames(reshape2::melt(LAA_data_1), c('year', 'age', 'LAA'))
data_1$type = 'Linf iid_y'
# Model 2:
LAA_data_2 = my_model_2$rep$LAA
rownames(LAA_data_2) = min_year:max_year
colnames(LAA_data_2) = 1:n_ages
data_2 = setNames(reshape2::melt(LAA_data_2), c('year', 'age', 'LAA'))
data_2$type = 'LAA iid'
# Model 3:
LAA_data_3 = my_model_3$rep$LAA
rownames(LAA_data_3) = min_year:max_year
colnames(LAA_data_3) = 1:n_ages
data_3 = setNames(reshape2::melt(LAA_data_3), c('year', 'age', 'LAA'))
data_3$type = 'Ecov Linf'
#Merge data:
plot_data = rbind(data_0, data_1, data_2, data_3)
# Plot:
ggplot(data = plot_data, aes(x = year, y = LAA, color = factor(type))) +
  geom_line() +
  theme_bw() +
  xlab(NULL) +
  ylab('Mean length-at-age') +
  labs(color = 'Model') +
  theme(legend.position = "bottom", legend.direction="horizontal") +
  facet_wrap(. ~ factor(age), nrow = 4, scales = 'free_y')
ggsave(filename = 'compare_LAA_Ex1.png', width = 190, height = 120, units = 'mm', dpi = 500)


## -----------------------------------------------------------------------------------------------------------------------
mycols = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

data0 = data.frame(model = 'OM', year = min_year:max_year, 
                   est = SS_report$timeseries$SpawnBio[SS_report$timeseries$Yr >= min_year &
                                                         SS_report$timeseries$Yr <= max_year],
                   sd = 0) # true
# Model 1:
this_model = my_model_1
model_name = 'Linf iid_y'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
data1 = cbind(model = model_name, year=26:75,
              dplyr::filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
# Model 2:
this_model = my_model_2
model_name = 'LAA iid'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
data2 = cbind(model = model_name, year=26:75,
              dplyr::filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
# Model 3:
this_model = my_model_3
model_name = 'Ecov Linf'
tmp = data.frame(name = names(this_model$sdrep$value),
                 est = this_model$sdrep$value, sd = this_model$sdrep$sd)
data3 = cbind(model = model_name, year=26:75,
              dplyr::filter(tmp, name=='log_SSB') %>% dplyr::select(-name))
# Merge data:
plot_data = rbind(data1, data2, data3)

# Make plot:
ggplot(plot_data, aes(year, exp(est), ymin=exp(est-1.96*sd), ymax=exp(est+1.96*sd),
                      fill=model, color=model)) +
  ylim(0,NA) + labs(y='SSB') +
  geom_ribbon(alpha=.3, color = NA) + geom_line(lwd=1) +
  geom_line(data = data0, aes(x = year, y = est)) +
  labs( color=NULL, fill=NULL) +
  scale_fill_manual(values = mycols[c(2,4,6,1)]) +
  scale_color_manual(values = mycols[c(2,4,6,1)]) +
  theme_bw() +
  theme(legend.position='top') 


## ---- eval=F------------------------------------------------------------------------------------------------------------
## mod$rep$logit_selpars # mean sel pars
## mod$rep$sel_repars # if time-varying selectivity turned on
## mod$rep$selAA # selectivity-at-age by block
## mod$sdrep # look for sel pars with NaN standard errors


## ---- eval=F------------------------------------------------------------------------------------------------------------
## 		# fix selectivity at 1 for ages 4-5 / 4 / 2-4 in blocks 1-3
## 		input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2,
## 					selectivity=list(model=rep("age-specific",3), re=sel_re[[m]],
## 						initial_pars=list(c(0.1,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)),
## 						fix_pars=list(4:5,4,2:4)),
## 					NAA_re = list(sigma='rec+1',cor='iid'))


## ---- eval=F------------------------------------------------------------------------------------------------------------
## mod.list <- file.path(getwd(),paste0("m",1:n.mods,".rds"))
## mods <- lapply(mod.list, readRDS)
## sapply(mods, function(x) check_convergence(x))
## sapply(mods, function(x) x$opt$obj) # get NLL
## lapply(mods, function(x) x$rep$logit_selpars)
## lapply(mods, function(x) x$rep$sel_repars)

