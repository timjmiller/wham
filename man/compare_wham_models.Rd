% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_wham_models.R
\name{compare_wham_models}
\alias{compare_wham_models}
\title{Compare multiple WHAM (or ASAP) models}
\usage{
compare_wham_models(
  mods,
  do.table = TRUE,
  do.plot = TRUE,
  fdir = getwd(),
  compare.opts = NULL,
  table.opts = NULL,
  plot.opts = NULL,
  fname = NULL,
  sort = NULL,
  calc.rho = NULL,
  calc.aic = NULL,
  do.print = NULL
)
}
\arguments{
\item{mods}{(named) list of fit WHAM/ASAP models. To read in ASAP model output, use \code{\link{read_asap3_fit}}. If no names are given, 'm1', 'm2', ...
will be used.}

\item{do.table}{T/F, produce table of AIC and/or Mohn's rho? Default = TRUE.}

\item{do.plot}{T/F, produce plots? Default = TRUE.}

\item{fdir}{character, path to directory to save table and/or plots. Default = \code{getwd()}.}

\item{compare.opts}{list of options to generate comparison results:
\describe{
  \item{\code{$stock}}{integer, which stock to include in results. Default = 1.}
  \item{\code{$region}}{integer, which region to include in results. Default = 1.}
}}

\item{table.opts}{list of options for AIC/rho table:
\describe{
  \item{\code{$fname}}{character, filename to save CSV results table (.csv will be appended). Default = \code{'model_comparison'}.}
  \item{\code{$sort}}{T/F, sort by AIC? Default = TRUE.}
  \item{\code{$calc.rho}}{T/F, calculate Mohn's rho? Retrospective analysis must have been run for all modes. Default = TRUE.}
  \item{\code{$calc.aic}}{T/F, calculate AIC? Default = TRUE.}
  \item{\code{$print}}{T/F, print table to console? Default = TRUE.}
  \item{\code{$save.csv}}{T/F, save table as a CSV file? Default = FALSE.}
}}

\item{plot.opts}{list of options for plots:
\describe{
  \item{\code{$out.type}}{character, either \code{'pdf'} or \code{'png'} (default = \code{'png'} because I am not sure \code{system('pdftk')} will work across platforms.)}
  \item{\code{$ci}}{vector of T/F, length = 1 (applied to all models) or number of models}
  \item{\code{$years}}{vector, which years to plot? Default = all (model and projection years).}
  \item{\code{$which}}{vector, which plots to make? Default = all. See details.}
  \item{\code{$relative.to}}{character, name of "base" model to plot differences relative to.}
  \item{\code{$alpha}}{scalar, (1-alpha)\% confidence intervals will be plotted. Default = 0.05 for 95\% CI.}
  \item{\code{$ages.lab}}{vector, overwrite model age labels.}
  \item{\code{$kobe.yr}}{integer vector (length = 1 (applied to all models or number of models, which year to use in Kobe plot (relative status). Default = terminal model year(s).}
  \item{\code{$M.age}}{integer, which age to use in M time-series plot. Default = \code{max(data$which_F_age)} (max age of F to use for full total F).}
  \item{\code{$return.ggplot}}{T/F, return a list of ggplot2 objects for later modification? Default = TRUE.}
  \item{\code{$kobe.prob}}{T/F, print probabilities for each model in each quadrant of Kobe plot? Default = TRUE.}
  \item{\code{$refpt}}{"XSPR" or "MSY", which reference point to use. Default = "XSPR".}
  \item{\code{$browse}}{Open html document in web browser (if $out.type = \code{'png'}) (default = TRUE).}
}}
}
\value{
a list with the following components:
  \describe{
    \item{\code{daic}}{Vector of delta-AIC by model (if \code{do.table=T} and \code{table.opts$calc.aic=T})}
    \item{\code{aic}}{Vector of AIC by model (if \code{do.table=T} and \code{table.opts$calc.aic=T})}
    \item{\code{rho}}{Matrix of Mohn's rho by model (if \code{do.table=T} and \code{table.opts$calc.rho=T})}
    \item{\code{best}}{Name of best model (lowest AIC) (if \code{do.table=T} and \code{table.opts$calc.aic=T})}
    \item{\code{tab}}{Results table of AIC and Mohn's rho (if \code{do.table=T})}
    \item{\code{g}}{List of ggplot2 objects for later modification (if \code{do.plot=T} and \code{plot.opts$return.ggplot=T})}
  }
}
\description{
After fitting multiple WHAM (or ASAP) models, \code{compare_wham_models} produces plots
and a table of AIC and Mohn's rho to aid model comparison.
}
\details{
\code{plot.opts$which} specifies which plots to make:
  \describe{
    \item{1}{3-panel of SSB (spawning stock biomass), F (fully-selected fishing mortality), and Recruitment}
    \item{2}{CV (coefficient of variation) for SSB, F, and Recruitment}
    \item{3}{Fleet selectivity (by block, averaged across years)}
    \item{4}{Index selectivity (by block, averaged across years)}
    \item{5}{Selectivity tile (fleets + indices, useful for time-varying random effects)}
    \item{6}{M time series (natural mortality, can specify which age with plot.opts$M.age)}
    \item{7}{M tile (useful for time-varying random effects)}
    \item{8}{3-panel of F X\% SPR, SSB at F_X\%SPR, and yield at F_X\%SPR}
    \item{9}{2-panel of relative status (SSB / SSB at F_X\%SPR and F / F_X\%SPR)}
    \item{10}{Kobe status (relative SSB vs. relative F)}
  }
If \code{plot.opts$return.ggplot = TRUE}, a list \code{g} is returned holding the above ggplot2 objects for later modification.
\code{g[[i]]} holds the plot corresponding to \code{i} above, e.g. \code{g[[2]]} is the CV plot.
}
\examples{
\dontrun{
base <- read_asap3_fit()
m1 <- fit_wham(input1)
m2 <- fit_wham(input2)
mods <- list(base=base, m1=m1, m2=m2)
res <- compare_wham_models(mods)
}

}
\seealso{
\code{\link{fit_wham}}, \code{\link{read_asap3_fit}, \code{\link{read_wham_fit}}}
}
