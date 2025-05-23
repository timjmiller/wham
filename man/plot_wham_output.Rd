% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_wham_output.R
\name{plot_wham_output}
\alias{plot_wham_output}
\title{Plot WHAM output}
\usage{
plot_wham_output(
  mod,
  dir.main = getwd(),
  out.type = "html",
  res = 72,
  plot.opts = NULL
)
}
\arguments{
\item{mod}{output from \code{\link{fit_wham}}}

\item{dir.main}{character, directory to save plots to (default = \code{getwd()})}

\item{out.type}{character, either \code{'html'}, \code{'pdf'}, or \code{'png'} (default = \code{'html'})}

\item{res}{resolution to save .png files (dpi)}

\item{plot.opts}{(optional) list of plot modifications}
}
\description{
Generates many output plots and tables for a fit WHAM model.
}
\details{
\code{out.type = 'html'} (default) makes a html file for viewing plot .png files and html tables of parameter estimates in a browser.
\code{out.type = 'pdf'} makes one pdf file of all plots and tables. 
\code{out.type = 'png'} creates a subdirectory `plots_png`` in \code{dir.main} and saves .png files within.
\code{out.type = 'pdf' or 'png'} makes LaTeX and pdf files of tables of parameter estimates.
(tabs: 'input data', 'diagnostics', 'results', 'ref_points', 'retro', and 'misc').

\code{plot.opts} holds optional arguments to modify plots:
  \describe{
    \item{\code{$ages.lab}}{Character vector, will change age labels in plots (default is \code{1:n.ages}).}
    \item{\code{$font.family}}{Font family, e.g. \code{"Times"}.}
    \item{\code{$browse}}{T/F whether to open the html file in a browser. Default = T.}
  }

Plot functions are located in \code{wham_plots_tables.R}
Table function is located in \code{par_tables_fn.R}
}
\examples{
\dontrun{
data("input4_SNEMAYT") # load fit wham model
mod <- fit_wham(input4_SNEMAYT)
plot_wham_output(mod)
}
}
\seealso{
\code{\link{fit_wham}}, \code{\link{wham_html}}, \code{wham_plots_tables}
}
