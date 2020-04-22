# WHAM example 6: Numbers-at-age options

# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref='naa')
# devtools::load_all()
# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex6_NAA.R"); etime <- Sys.time(); runtime = etime - btime;
# 7.6 min

context("Ex 6: Numbers-at-age")

test_that("Ex 6 works",{
path_to_examples <- system.file("extdata", package="wham")
ex5_test_results <- readRDS(file.path(path_to_examples,"ex5_test_results.rds"))

# these 4 models are in Haikun's paper
# NAA_re = list(cor='iid', sigma='rec+1')
# NAA_re = list(cor='ar1_a', sigma='rec+1')
# NAA_re = list(cor='ar1_y', sigma='rec+1')
# NAA_re = list(cor='2dar1', sigma='rec+1')

# NAA_re = list(cor='iid', sigma='rec') # previously "random_recruitment"
# NAA_re = list(cor='ar1_a', sigma='rec') # new option
# NAA_re = list(cor='ar1_y', sigma='rec') # new option
# NAA_re = list(cor='2dar1', sigma='rec') # not allowed (error message)

NAA_re = NULL # SCAA, default

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
df.mods <- data.frame(NAA_cor = c('---','iid','ar1_a','ar1_y','2dar1','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('---',rep("rec",4),rep("rec+1",4)), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
# df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

mods <- vector("list",n.mods)
for(m in 1:n.mods){
  # set up environmental covariate data and model options
  # see ?prepare_wham_input
  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    # lag = 1, # GSI in year t affects Rec in year t + 1
    lag = 0, # GSI in year t affects M in same year
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    where = "M", # GSI affects natural mortality
    how = ifelse(df.mods$Ecov_link[m]==0,0,1), # 0 = no effect (but still fit Ecov to compare AIC), 1 = mean
    link_model = c(NA,"linear","poly-2")[df.mods$Ecov_link[m]+1])

  m_model <- df.mods$M_model[m]
  if(df.mods$M_model[m] == '---') m_model = "age-specific"
  if(df.mods$M_model[m] %in% c("constant","weight-at-age")) est_ages = 1
  if(df.mods$M_model[m] == "age-specific") est_ages = 1:asap3$dat$n_ages
  if(df.mods$M_model[m] == '---') est_ages = NULL
  M <- list(
    model = m_model,
    re = df.mods$M_re[m],
    est_ages = est_ages
  )
  if(m_model %in% c("constant","weight-at-age")) M$initial_means = 0.28

  # Generate wham input from ASAP3 and Ecov data
  input <- prepare_wham_input(asap3, recruit_model = 2,
                              model_name = "Ex 5: Yellowtail Flounder with GSI effects on M",
                              ecov = ecov,
                              selectivity=list(model=rep("logistic",6),
                                               initial_pars=c(rep(list(c(3,3)),4), list(c(1.5,0.1), c(1.5,0.1))),
                                               fix_pars=c(rep(list(NULL),4), list(1:2, 1:2))),
                              M=M)

  # overwrite age comp model (all models use logistic normal)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # Full state-space model, abundance is the state vector
  input$data$use_NAA_re = 1
  input$data$random_recruitment = 0
  input$map = input$map[!(names(input$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))] # remove from map (i.e. estimate)
  input$map$log_R = factor(rep(NA, length(input$par$log_R)))
  input$random = c(input$random, "log_NAA")

  # Fit model
  mod <- suppressWarnings(fit_wham(input, do.retro=F, do.osa=F))

  expect_equal(as.numeric(mod$opt$par), ex5_test_results$pars[[m]], tolerance=1e-3) # parameter values
  expect_equal(as.numeric(mod$opt$obj), ex5_test_results$nll[m], tolerance=1e-6) # nll
}

})
