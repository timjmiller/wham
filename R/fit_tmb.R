#' Fit TMB model using nlminb
#'
#' Runs optimization on the TMB model using \code{\link[stats:nlminb]{stats::nlminb}}.
#' If specified, takes additional Newton steps and calculates standard deviations.
#' Internal function called by \code{\link{fit_wham}}.
#'
#' @param model Output from \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}.
#' @param n.newton Integer, number of additional Newton steps after optimization. Default = \code{3}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#' @param do.check T/F, check if model parameters are identifiable? Runs \code{\link[TMBhelper::Check_Identifiable]{TMBhelper::Check_Identifiable}}. Default = \code{TRUE}.
#' @return \code{model}, appends the following:
#'   \describe{
#'     \item{\code{model$opt}}{Output from \code{\link[stats:nlminb]{stats::nlminb}}}
#'     \item{\code{model$date}}{System date}
#'     \item{\code{model$dir}}{Current working directory}
#'     \item{\code{model$rep}}{model$report()}
#'     \item{\code{model$TMB_version}}{Version of TMB installed}
#'     \item{\code{model$parList}}{List of parameters, \code{model$env$parList()}}
#'     \item{\code{model$final_gradient}}{Final gradient, \code{model$gr()}}
#'     \item{\code{model$sdrep}}{Estimated standard deviations for model parameters, \code{\link[TMB:sdreport]{TMB::sdreport}}}
#'   }
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}, \code{\link{TMBhelper::Check_Identifiable}}
#'
fit_tmb = function(model, n.newton=3, do.sdrep=TRUE, do.check=TRUE)
{
  model$opt <- stats::nlminb(model$par, model$fn, model$gr, control = list(iter.max = 1000, eval.max = 1000))

  if(n.newton){ # Take a few extra newton steps
    tryCatch(for(i in 1:n.newton) { 
      g <- as.numeric(model$gr(model$opt$par))
      h <- stats::optimHess(model$opt$par, model$fn, model$gr)
      model$opt$par <- model$opt$par - solve(h, g)
      model$opt$objective <- model$fn(model$opt$par)
    }, error = function(e) {err <<- conditionMessage(e)}) # still want fit_tmb to return model if newton steps error out
  }
  if(exists("err")) model$err <- err # store error message to print out in fit_wham

  if(do.check){
    ParHat = model$env$last.par.best[-c(model$env$random)]
    Gr = model$gr(ParHat)
    if(any(Gr > 0.01)){
      df <- data.frame(param = names(ParHat),
                       MLE = ParHat,
                       gr.at.MLE = Gr)
      ind.hi <- which(Gr > 0.01)
      model$badpar <- df[ind.hi,]
      warning(paste("","Some parameter(s) have high gradients at the MLE:","",
        paste(capture.output(print(model$badpar)), collapse = "\n"), sep="\n"))
    } else {
      test <- TMBhelper::Check_Identifiable(model)
      if(length(test$WhichBad) > 0){
        bad.par <- as.character(test$BadParams$Param[test$BadParams$Param_check=='Bad'])
        bad.par.grep <- grep(bad.par, test$BadParams$Param)
        model$badpar <- test$BadParams[bad.par.grep,]
        warning(paste("","Some fixed effect parameter(s) are not identifiable, consider removing",
          "them from the model by setting input$par at their MLE and input$map = NA.","",
          paste(capture.output(print(test$BadParams[bad.par.grep,])), collapse = "\n"), sep="\n"))    
      }
    }
  }

  model$date = Sys.time()
  model$dir = getwd()
  model$rep <- model$report()
  model$TMB_version = packageVersion("TMB")
  model$parList = model$env$parList()
  model$final_gradient = model$gr()

  # only do sdrep if no error
  if(do.sdrep & !exists("err"))
  {
    model$sdrep <- try(TMB::sdreport(model))
    model$is_sdrep = !is.character(model$sdrep)
    if(model$is_sdrep) model$na_sdrep = any(is.na(summary(model$sdrep,"fixed")[,2]))
    else model$na_sdrep = NA
  } else {
    model$is_sdrep = FALSE
    model$na_sdrep = NA
  }

  return(model)
}
