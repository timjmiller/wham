# WHAM example 5: Ecov and age-year effects on natural mortality

# devtools::load_all()
# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex5_M.R"); etime <- Sys.time(); runtime = etime - btime;
# 7.6 min

context("Ex 5: Natural mortality")

test_that("Ex 5 works",{
path_to_examples <- system.file("extdata", package="wham")
ex5_test_results <- readRDS(file.path(path_to_examples,"ex5_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)

df.mods <- data.frame(M_model = c(rep("---",3),"age-specific","weight-at-age",rep("constant",6),"age-specific","age-specific",rep("constant",3),"---"),
                      M_re = c(rep("none",6),"ar1_y","2dar1","none","none","2dar1","none","2dar1",rep("ar1_a",3),"2dar1"),
                      Ecov_process = rep("ar1",17),
                      Ecov_link = c(0,1,2,rep(0,5),1,2,1,2,2,0,1,2,0), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)

mods <- vector("list",n.mods)
for(m in c(1:3,5:10)){
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