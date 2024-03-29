% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_peel.R
\name{fit_peel}
\alias{fit_peel}
\title{Fit model peeling off \emph{i} years of data}
\usage{
fit_peel(
  peel,
  input,
  do.sdrep = FALSE,
  n.newton = 3,
  MakeADFun.silent = FALSE,
  retro.silent = FALSE,
  save.input = FALSE
)
}
\arguments{
\item{peel}{Integer, number of years of data to remove before model fitting.}

\item{input}{input with same structure as that provided by \code{\link{prepare_wham_input}}. May want to use input$par = model$parList to start at MLEs.}

\item{do.sdrep}{T/F, calculate standard deviations of model parameters? Default = \code{FALSE}.}

\item{n.newton}{integer, number of additional Newton steps after optimizafit_tmbtion for each peel. Default = \code{3}.}

\item{MakeADFun.silent}{T/F, Passed to silent argument of \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}. Default = \code{FALSE}.}

\item{retro.silent}{T/F, Passed to argument of internal fit_peel function. Determines whether peel number is printed to screen. Default = \code{FALSE}.}

\item{save.input}{T/F, should modified input list be saved? Necessary to project from a peel but increases model object size. Default = \code{FALSE}.}
}
\value{
\code{out}, output of \code{\link{fit_tmb}} for peel \emph{i}
}
\description{
Internal function called by \code{\link{retro}} for \emph{i} in 1--\code{n.peels}.
Fits the model peeling off \emph{i} years of data (calls \code{\link{fit_tmb}}).
}
\seealso{
\code{\link{fit_wham}}, \code{\link{retro}}, \code{\link{fit_tmb}}
}
