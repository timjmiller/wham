% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_M.R
\name{set_M}
\alias{set_M}
\title{Specify model and parameter configuration for natural mortality}
\usage{
set_M(input, M)
}
\arguments{
\item{input}{list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})}

\item{M}{(optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see details)

\code{M} specifies estimation options for natural mortality and can overwrite M-at-age values specified in the ASAP data file.
If \code{NULL}, the M-at-age matrix from the ASAP data file is used (M fixed, not estimated). To estimate M-at-age
shared/mirrored among some but not all ages, modify \code{M$means_map} (see vignette for more details). \code{M} is a list 
with the following entries:
  \describe{
    \item{$mean_model}{Character describing the type of model for M stock and regional models for natural mortality. Options are:
      \describe{
        \item{"fixed-M"}{Use initial values from ASAP3 dat files or \code{$initial_means} for (mean) M as fixed values. If no ASAP3 files
          and \code{$initial_means} is not provided, default is M = 0.2 for all stocks, regions and ages}
        \item{"estimate-M"}{estimate one or more (mean) M parameters. Default is to estimate a single M shared across all stocks and ages, but
          use \code{$means_map} to fix or estimate parameters for specific stocks, regions, ages.}
        \item{"weight-at-age"}{specifies M as a function of weight-at-age, \eqn{M_y,a = exp(b0 + b1*log(W_y,a))}, as in
          \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)} and
          \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035}{Miller & Hyun (2018)}.
          Default is to estimate a single model shared across all stocks and regions, but
          use \code{$means_map[s,r,1]} to fix or estimate the intercept for specific stocks, regions. See also \code{$logb_prior}
          and \code{$initial_b} configuring the slope on log scale.}
      }
    }
    \item{$initial_means}{array (n_stocks x n_regions x n_ages) of initial/mean M by stock, region and age. If \code{NULL}, initial 
      mean M-at-age values for a given stock and region are taken from the first row of the MAA matrix in the ASAP data file. If no
       ASAP data file, M = 0.2 is the default. If \code{$mean_model} is "weight-at-age" only 
      elements for the first age (\code{$initial_means[,,1]}) are used (for the intercept of log(M)).}
    \item{$means_map}{array (n_stocks x n_regions x n_ages) of NA or integers ( 0 <= max <= n_stocks * n_regions * n_ages) indicating 
      which ages to estimate (mean) M and whether to set any ages to be identical. E.g. in a model with 2 stock, 2 regions 
      and 6 ages with constant M estimated for each stock across regions and ages  \code{$M_ages_map[1,,] = 1} 
      and \code{$M_ages_map[2,,] = 2}. \code{$M_ages_map[1,1,] = c(NA,1,1,2,2,3)} will fix M for age 1 at the initial value, 
      and estimates for ages 2 and 3 are identical as are those for ages 4 and 5 and different from age 6+ for stock 1 and 
      region 1. If \code{NULL}, specifies all ages fixed at \code{M$initial_means}.  If \code{$mean_model} is "weight-at-age"
      these are used for all stocks and regions and only the elements for the first age (\code{$M_ages_map[,,1]}) 
      are used (for the intercept of log(M)).}
    \item{$intial_MAA}{array (n_stocks x n_regions x n_years x n_ages) of initial values for M at age. Intended to be uses when nothing pertaining to M estimated.}
    \item{$b_model}{"constant","stock","region", "stock-region" defining whether parameter is constant, stock-specific, region-specific, 
      stock- and region-specific. Only used if \code{M$mean_model} = "weight-at-age".}
    \item{$b_prior}{T/F, should a N(mu, 0.08) prior (where mu = log(0.305) by default) be used on log_b? Based on Fig. 1 and Table 1 
      (marine fish) in  \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)}. (Only used if 
      \code{$mean_model} is "weight-at-age").}
    \item{$intial_b}{if any elements of \code{$mean_model} is "weight-at-age", initial value for mean b for weight-at-age
      model.}
    \item{$re_model}{Character matrix (n_stocks x n_regions) of options for time- and age-varying (random effects) on M by stock and region.
      Possible values are:
      \describe{
        \item{"none"}{(default) No random effects by age or year.}
        \item{"iid_a"}{uncorrelated M by age, constant in time.}
        \item{"iid_y"}{uncorrelated M by year, constant all ages.}
        \item{"ar1_a"}{M correlated by age (AR1), constant in time.}
        \item{"ar1_y"}{M correlated by year (AR1), constant all ages.}
        \item{"iid_ay"}{M uncorrelated by year and age (2D).}
        \item{"ar1_ay"}{M correlated by year and age (2D AR1), as in \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047}{Cadigan
          (2016)}.}
      }
    }
    \item{$re_map}{array (n_stocks x n_regions x n_ages) of NA and integers (1 <= max <= n_ages) indicating which ages, for a given 
      stock and region, have random effects (not NA) and whether to set RE for any ages to be identical. E.g. in a model with 2 stock, 
      2 regions and 6 ages, \code{$re_map[2,1,] = c(NA,1,1,2,2,3)} will not estimate RE for age 1, and those for ages 
      2 and 3 are identical as are those for ages 4 and 5 and different from age 6+ for stock 2 and region 1. If \code{NULL}, 
      and \code{$re_model} specifies M random effects at age, at least two ages must be 
      specified for correlation among ages to be estimated.}
    \item{$sigma_vals}{n_stocks x n_regions matrix Initial standard deviation value to use for the M random effects. Values are not used 
      if \code{M$re_model} = "none". Otherwise, a single value. If unspecified all values are 0.1.}
    \item{$cor_vals}{n_stocks x n_regions x 2 array of initial correlation values to use for the M deviations. If unspecified all initial 
      values are 0. When \code{M$re_model} = 
      \describe{
        \item{"iid_a", "iid_y", "iid_ay" or "none"}{values are not used.}
        \item{"ar1_a" }{first value cor_vals[s,r,1] is used.}
        \item{"ar1_y" }{second value cor_vals[s,r,2] is used.}
        \item{"ar1_ay"}{First is for "age", second is for "year".}
      }
    }
    \item{$sigma_map}{n_stocks x n_region matrix of NA or integers indicating which random effects sd is estimated and whether to 
      set any to be identical. If not supplied a single sd will be estimated for any stock and region where $re_model is other than "none".}
    \item{$cor_map}{n_stocks x n_region matrix x 2 array of NA or integers indicating which random effects correlation parameters are estimated
      and whether to set any to be identical. If not supplied a single value for age and/or year will be estimated for any stock and region where 
      $re_model is other than "none", "iid_a", "iid_y".}
  }}
}
\value{
a named list with same elements as the input provided with natural mortality options modified.
}
\description{
Specify model and parameter configuration for natural mortality
}
\examples{
\dontrun{
wham.dir <- find.package("wham")
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3)
M = list(mean_model = "estimate-M")
input <- set_q(input, M = M) #estimate a constant M parameters
}

}
\seealso{
\code{\link{prepare_wham_input}}
}
