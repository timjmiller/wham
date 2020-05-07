# WHAM example 6: Numbers-at-age options

# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref='naa')
# library(wham)
# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex6_NAA.R"); etime <- Sys.time(); runtime = etime - btime;
# 14 min

context("Ex 6: Numbers-at-age")

test_that("Ex 6 works",{
path_to_examples <- system.file("extdata", package="wham")
ex6_test_results <- readRDS(file.path(path_to_examples,"ex6_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)

df.mods <- data.frame(NAA_cor = c('---','iid','ar1_y','iid','ar1_a','ar1_y','2dar1','iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('---',rep("rec",2),rep("rec+1",4),rep("rec",2),rep("rec+1",4)),
                      GSI_how = c(rep(0,7),rep(2,6)), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
# df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

fit.mods <- c(1:4,6:13)
for(m in fit.mods){
  NAA_list <- list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"])
  if(NAA_list$sigma == '---') NAA_list = NULL

  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    lag = 1, # GSI in year t affects Rec in year t + 1
    process_model = 'ar1', # "rw" or "ar1"
    where = "recruit", # GSI affects recruitment
    how = df.mods$GSI_how[m], # 0 = no effect (but still fit Ecov to compare AIC), 2 = limiting
    link_model = "linear")

  recruit = ifelse(m == 1,2,3)
  input <- prepare_wham_input(asap3, recruit_model = recruit, # Bev Holt recruitment
                              model_name = "Ex 6: Numbers-at-age",
                              selectivity=list(model=rep("age-specific",3), 
                                re=rep("none",3), 
                                initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,0.5,0.5,0.5,0.5)), 
                                fix_pars=list(c(4,5),4,2)),
                              NAA_re = NAA_list,
                              ecov=ecov)

  # overwrite age comp model (all models use logistic normal)
  input$data$age_comp_model_indices = rep(7, input$data$n_indices)
  input$data$age_comp_model_fleets = rep(7, input$data$n_fleets)
  input$data$n_age_comp_pars_indices = rep(1, input$data$n_indices)
  input$data$n_age_comp_pars_fleets = rep(1, input$data$n_fleets)
  input$par$index_paa_pars = rep(0, input$data$n_indices)
  input$par$catch_paa_pars = rep(0, input$data$n_fleets)
  input$map = input$map[!(names(input$map) %in% c("index_paa_pars", "catch_paa_pars"))]

  # Fit model
  mod <- suppressWarnings(fit_wham(input, do.retro=F, do.osa=F))

  cat("---------------------------------------------------------------------------------------------------------------------\n")
  cat(paste0("Model ",m,"\n"))
  cat("---------------------------------------------------------------------------------------------------------------------\n")
  # expect_equal(as.numeric(mod$opt$par), ex6_test_results$pars[[m]], tolerance=1e-3) # parameter values
  expect_equal(as.numeric(mod$opt$obj), ex6_test_results$nll[m], tolerance=1e-3) # nll
}

})