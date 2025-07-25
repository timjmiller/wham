% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mohns_rho.R
\name{mohns_rho}
\alias{mohns_rho}
\title{Calculate Mohn's rho for a WHAM model with peels}
\usage{
mohns_rho(model)
}
\arguments{
\item{model}{A fit WHAM model, output from \code{\link{fit_wham}} with \code{do.retro = TRUE}.}
}
\value{
\code{rho}, a vector of Mohn's rho
}
\description{
Calculate Mohn's rho for a WHAM model with peels
}
\examples{
\dontrun{
data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
mod <- fit_wham(input4_SNEMAYT) # using default values: do.retro = T, n.peels = 7
mohns_rho(mod) # calculate Mohn's rho
}
}
\seealso{
\code{\link{fit_wham}}, \code{\link{retro}}
}
