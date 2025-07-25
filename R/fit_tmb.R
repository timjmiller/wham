#' Fit TMB model using nlminb
#'
#' Runs optimization on the TMB model using \code{\link[stats:nlminb]{stats::nlminb}}.
#' If specified, takes additional Newton steps and calculates standard deviations.
#' Internal function called by \code{\link{fit_wham}}.
#'
#' @param model Output from \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}.
#' @param n.newton Integer, number of additional Newton steps after optimization. Default = \code{3}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{TMB::sdreport}}. Default = \code{TRUE}.
#' @param do.check T/F, check if model parameters are identifiable? Runs internal \code{\link{check_estimability}}, originally provided by the [TMBhelper](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper) package. Default = \code{TRUE}.
#' @param save.sdrep T/F, save the full \code{\link[TMB]{TMB::sdreport}} object? If \code{FALSE}, only save \code{\link[TMB:summary.sdreport]{summary.sdreport)}} to reduce model object file size. Default = \code{FALSE}.
#' @param use.optim T/F, use \code{\link[stats]{stats::optim}} instead of \code{\link[stats]{stats::nlminb}}? Default = \code{FALSE}.
#' @param opt.control list of control parameters to pass to optimizer. For nlminb default = list(iter.max = 1000, eval.max = 1000). For optim default = list(maxit=1000).
#' @return \code{model}, appends the following:
#'   \describe{
#'     \item{\code{model$opt}}{Output from \code{\link[stats:nlminb]{stats::nlminb}}}
#'     \item{\code{model$date}}{System date}
#'     \item{\code{model$dir}}{Current working directory}
#'     \item{\code{model$rep}}{model$report(model$env$last.par.best)}
#'     \item{\code{model$TMB_version}}{Version of TMB installed}
#'     \item{\code{model$parList}}{List of parameters, \code{model$env$parList(x = model$opt$par, par = model$env$last.par.best)}}
#'     \item{\code{model$final_gradient}}{Final gradient, \code{model$gr(model$opt$par)}}
#'     \item{\code{model$sdrep}}{Estimated standard deviations for model parameters, \code{\link[TMB:sdreport]{TMB::sdreport}} or \code{\link[TMB:summary.sdreport]{summary.sdreport)}}}
#'   }
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}, \code{check_estimability}
#'
#'
#' @export
#' @md
fit_tmb <- function(model, n.newton=3, do.sdrep=TRUE, do.check=FALSE, save.sdrep=FALSE, use.optim=FALSE, opt.control = NULL)
{
  if(use.optim){
    if(is.null(opt.control)) opt.control <- list(maxit = 1000)
    message("Using stats::optim for optimization rather than stats::nlminb with these control parameters:")
    message(paste0(names(opt.control), ": ", opt.control, "\n"))
    tryCatch(model$opt <- stats::optim(model$par, model$fn, model$gr, control = opt.control), 
    error = function(e) {model$opt_err <<- conditionMessage(e)})
    model$opt$objective <- model$opt$value #try to make sure calls to get nll will work for both optimizers
  }
  else {
    if(is.null(opt.control)) opt.control <- list(iter.max = 1000, eval.max = 1000)
    message("Using stats::nlminb for optimization with these control parameters:")
    message(paste0(names(opt.control), ": ", opt.control, "\n"))
    tryCatch(model$opt <- stats::nlminb(model$par, model$fn, model$gr, control = opt.control), 
    error = function(e) {model$opt_err <<- conditionMessage(e)})
  }
  if(is.null(model$opt_err)){
    Gr <- model$gr(model$opt$par)
    opt_before <- model$opt
    if(n.newton) if(!any(is.na(Gr))) if(max(abs(Gr))<1){ # Take a few extra newton steps when useful
      tryCatch(for(i in 1:n.newton) { 
        g <- as.numeric(model$gr(model$opt$par))
        h <- stats::optimHess(model$opt$par, model$fn, model$gr)
        model$opt$par <- model$opt$par - solve(h, g)
        model$opt$objective <- model$fn(model$opt$par)
      }, error = function(e) {model$err <<- conditionMessage(e)}) # still want fit_tmb to return model if newton steps error out
      model$newton_info <- list()
      model$newton_info$nll_diff <- model$opt$obj - opt_before$obj
      model$newton_info$par_diff <- model$opt$par - opt_before$par
      model$newton_info$grad_diff <- model$gr(model$opt$par) - Gr
      if(model$newton_info$nll_diff > 0) message(paste0("Note: Newton steps resulted in an increase in the NLL of ", model$newton_info$nll_diff,"."))
    }
    is.re <- length(model$env$random)>0
    fe <- model$opt$par
    if(is.re) model$env$last.par.best[-c(model$env$random)] <- fe
    else  model$env$last.par.best <- fe
    
    Gr <- model$gr(fe)
    model$parList <- model$env$parList(x = fe, par = model$env$last.par.best)
    model$final_gradient <- Gr
    if(do.check){
      if(any(abs(Gr) > 0.01)){
        df <- data.frame(param = names(fe),
                         MLE = fe,
                         gr.at.MLE = Gr)
        ind.hi <- which(abs(Gr) > 0.01)
        model$badpar <- df[ind.hi,]
        warning(paste("","Some parameter(s) have high |gradients| at the MLE:","",
          paste(capture.output(print(model$badpar)), collapse = "\n"), sep="\n"))
      } else {
        test <- check_estimability(model)
        if(length(test$WhichBad) > 0){
          bad.par <- as.character(test$BadParams$Param[test$BadParams$Param_check=='Bad'])
          bad.par.grep <- grep(bad.par, test$BadParams$Param)
          model$badpar <- test$BadParams[bad.par.grep,]
          warning(paste("","Some fixed effect parameter(s) are not identifiable.",
            "Consider 1) removing them from the model by fixing input$par and input$map = NA, or",
            "2) changing your model configuration.","",
            paste(capture.output(print(test$BadParams[bad.par.grep,])), collapse = "\n"), sep="\n"))    
        }
      }
    }
  }


  model$date <- Sys.time()
  model$dir <- getwd()
  model$rep <- model$report(model$env$last.par.best)
  TMB_commit <- packageDescription("TMB")$GithubSHA1
  model$TMB_commit <- ifelse(is.null(TMB_commit), "local install", paste0("Github (kaskr/adcomp@", TMB_commit, ")")) 
  TMB_version <- packageDescription("TMB")$Version
  model$TMB_version <- paste0(TMB_version, " / ", model$TMB_commit, ")")


  if(do.sdrep) # only do sdrep if no error
  {
    model <- do_sdreport(model, save.sdrep = save.sdrep, TMB.bias.correct = FALSE, TMB.jointPrecision = FALSE)
  } else {
    model$is_sdrep <- FALSE
    model$na_sdrep <- NA
  }

  return(model)
}


#' Extract fixed effects
#' Originally provided by the [TMBhelper](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper) package.
#' Internal function called by \code{\link{check_estimability}}.
#'
#' \code{extract_fixed} extracts the best previous value of fixed effects, in a way that works for both mixed and fixed effect models
#'
#' @param obj, The compiled object
#'
#' @return A vector of fixed-effect estimates
#' @md
extract_fixed <- function( obj ){
  if( length(obj$env$random)==0 ){
    Return <- obj$env$last.par.best
  }else{
    Return <- obj$env$last.par.best[-c(obj$env$random)]
  }
  return( Return )
}


#' Check for identifiability of fixed effects
#' Originally provided by the [TMBhelper](https://github.com/kaskr/TMB_contrib_R/tree/master/TMBhelper) package.
#' Internal function called by \code{\link{fit_tmb}}.
#'
#' \code{check_estimability} calculates the matrix of second-derivatives of the marginal likelihood
#' w.r.t. fixed effects, to see if any linear combinations are not estimable (i.e. cannot be
#' uniquely estimated conditional upon model structure and available data, e.g., resulting
#' in a likelihood ridge and singular, non-invertable Hessian matrix)
#'
#' @param obj The compiled object
#' @param h optional argument containing pre-computed Hessian matrix
#'
#' @return A tagged list of the hessian and the message

#' @export
#' @md
check_estimability <- function( obj, h ){

  # Extract fixed effects
  ParHat <- extract_fixed( obj )

  # Check for problems
  Gr <- obj$gr( ParHat )
  if( any(abs(Gr)>0.01) ) stop("Some |gradients| are high, please improve optimization and only then use `check_estimability`")

  # Finite-different hessian
  List <- NULL
  if(missing(h)){
    List[["Hess"]] <- optimHess( par=ParHat, fn=obj$fn, gr=obj$gr )
  }else{
    List[["Hess"]] <- h
  }

  # Check eigendecomposition
  List[["Eigen"]] <- eigen( List[["Hess"]] )
  List[["WhichBad"]] <- which( List[["Eigen"]]$values < sqrt(.Machine$double.eps) )

  # Check result
  if( length(List[["WhichBad"]])==0 ){
    message( "All parameters are estimable" )
  }else{
    # Check for parameters
    RowMax <- apply( List[["Eigen"]]$vectors[,List[["WhichBad"]],drop=FALSE], MARGIN=1, FUN=function(vec){max(abs(vec))} )
    List[["BadParams"]] <- data.frame("Param"=names(obj$par), "MLE"=ParHat, "Param_check"=ifelse(RowMax>0.1, "Bad","OK"))
    message( List[["BadParams"]] )
  }

  return( invisible(List) )
}
