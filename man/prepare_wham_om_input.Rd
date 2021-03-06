% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_wham_om_input.R
\name{prepare_wham_om_input}
\alias{prepare_wham_om_input}
\title{Prepare input data and parameters for WHAM model from a few basic inputs to simulate populations and observations}
\usage{
prepare_wham_om_input(
  basic_info,
  model_name = "WHAM setup to simulate unnamed stock",
  recruit_model = 2,
  ecov = NULL,
  selectivity = NULL,
  M = NULL,
  NAA_re = NULL,
  proj.opts = NULL
)
}
\arguments{
\item{basic_info}{(optional) list basic inputs for data and parameters}

\item{model_name}{character, name of stock/model}

\item{recruit_model}{numeric, option to specify stock-recruit model (see details)}

\item{ecov}{(optional) named list of environmental covariate data and parameters (see details)}

\item{selectivity}{(optional) list specifying selectivity options by block: models, initial values, parameters to fix, and random effects (see details)}

\item{M}{(optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see details)}

\item{NAA_re}{(optional) list specifying random effect options on numbers-at-age (see details)}

\item{proj.opts}{(optional) list specifying projection options and used with \code{prepare_wham_om_proj} to modify input for projection years.}
}
\value{
a named list with the following components:
  \describe{
    \item{\code{data}}{Named list of data, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
    \item{\code{par}}{Named list of parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
    \item{\code{map}}{Named list defining how to optionally collect and fix parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
    \item{\code{random}}{Character vector of parameters to treat as random effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
    \item{\code{years}}{Numeric vector of years to fit WHAM model (specified in ASAP3 .dat file)}
    \item{\code{ages.lab}}{Character vector of age labels, ending with plus-group (specified in ASAP3 .dat file)}
    \item{\code{model_name}}{Character, name of stock/model (specified in call to \code{prepare_wham_input})}
  }
}
\description{
This function takes an optional \code{basic_info} list and prepares the data and parameter settings that can be used with
\code{\link[TMB:MakeADFun]{TMB::MakeADFun}} to generate a model to do operating model simulations.
By default, this will set up a SCAA version like \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP}.
}
\details{
\code{basic_info} may or may not be specified. If not, default values will be used to generate input and used in the function. 
It is also provided as output and may be modified and provided to the function again.  If it is specified, it has the following elements:
  \describe{
    \item{$modyears}{Years which the model (non-projection years) should span. Default is 1975:2014}
    \item{$maturity}{Vector of age-specific proportion mature. Default is 1/(1 + exp(-1*(1:na - 5))) where na is the number of age classes}
    \item{$fracyr_spawn}{fraction of the calendar year elapsed when spawning occurs}
    \item{$waa_catch}{Matrix of age-specfic weights of fish in catch by fleet (n_fleets x n_ages)}
    \item{$waa_totcatch}{Vector of age-specfic weights of fish in total catch across all fleets}
    \item{$waa_indices}{Matrix of age-specfic weights of fish observed by index (n_indices x n_ages)}
    \item{$waa_ssb}{Vector of age-specific weights of fish in population at the time of spawning}
    \item{$waa_jan1}{Vector of age-specific weights of fish in population at beginning of the calendar year}
    \item{$catch_cv}{Vector of fleet-specific CVs to use for aggregate catch}
    \item{$catch_Neff}{Vector of fleet-specific multinomial effective sample sizes for age compositions}
    \item{$index_cv}{Vector of index-specific CVs to use for aggregate indices}
    \item{$index_Neff}{Vector of index-specific multinomial effective sample sizes for age compositions}
    \item{$fracyr_indices}{Vector of fractions of the calendar year elapsed when the index observations are obtained}
    \item{$sel_model_fleets}{Vector of integers specifying selectivity models for each fleet. E.g. increasing logistic = 2}
    \item{$sel_model_indices}{Vector of integers specifying selectivity models for each index. E.g. increasing logistic = 2}
    \item{$q}{Vector of fully-selected catchabilities to use for each index}
    \item{$F}{Matrix of annual fully-selected fishing mortalities to use for each fleet (n_years x n_fleets)}
    \item{$M}{Vector of age-specific natural mortality rates to apply to the population. May be modified by specifying the \code{M} argument defined below.}
    \item{$mean_rec_pars}{Vector of recruitment parameters. Size depends on \code{recruit_model} argument defined below.}
    \item{$N1}{Vector of abundances at age for the population at the beginning of the first year of the model.}
    \item{$NAA_sigma}{Vector of process error standard deviations to use with numbers at age random effects}
    \item{$NAA_rho}{Vector of AR1 correlation parameters to use with numbers at age random effects. At most 2 values for correlation with age and/or year.}
    \item{$catch_paa_pars}{Vector of any dispersion parameters that may be used with alternative likelihood assumptions for fleet age composition observations.}
    \item{$index_paa_pars}{Vector of any dispersion parameters that may be used with alternative likelihood assumptions for index age composition observations.}
  }

\code{recruit_model} specifies the stock-recruit model. See \code{wham.cpp} for implementation.
  \describe{
    \item{= 1}{Random walk, i.e. predicted recruitment in year i = recruitment in year i-1}
    \item{= 2}{(default) Random about mean, i.e. steepness = 1}
    \item{= 3}{Beverton-Holt}
    \item{= 4}{Ricker}
  }
Note: by default a SCAA model is assumed (\code{NAA_re = NULL}), which parameterizes recruitment in every year as separate fixed effect parameters,
but in that case no stock-recruit function is allowed. A warning message is printed if \code{recruit_model > 2} and \code{NAA_re = NULL}.
If you wish to use a stock-recruit function for expected recruitment, choose recruitment deviations as random effects, 
either only age-1 (\code{NAA_re = list(sigma='rec')}) or all ages (\code{NAA_re = list(sigma='rec+1')}, "full state-space" model).
See below for details on \code{NAA_re} specification.

\code{ecov} specifies any environmental covariate data and model. Environmental covariate data need not span
the same years as the fisheries data. It can be \code{NULL} if no environmental data are to be fit.
Otherwise, it must be a named list with the following components:
  \describe{
    \item{$label}{Name(s) of the environmental covariate(s). Used in printing.}
    \item{$mean}{Mean observations (matrix). Missing values = NA.}
    \item{$logsigma}{Observation standard error (log). Options:
      \describe{
        \item{Matrix with same dimensions as \code{$mean}}{Specified values for each time step }
        \item{Single value per ecov, numeric vector or matrix w/ dim 1 x n.ecov}{Specified value shared among time steps}
        \item{Vector with estimate options for each ecov, length = n.ecov}{
          \code{'est_1'}: One value shared among time steps, estimated when used with \code{\link{fit_wham}}.
          \code{'est_re'}: Values for each time step are random effects with two parameters (mean, var)}
      }
    }
    \item{$year}{Years corresponding to observations (vector of same length as \code{$mean} and \code{$logsigma})}
    \item{$use_obs}{T/F (or 0/1) vector/matrix of the same dimension as \code{$mean} and \code{$logsigma}.
    Use the observation? Can be used to ignore subsets of the ecov without changing data files.}
    \item{$lag}{Offset between the ecov observations and their affect on the stock.
    I.e. if SST in year \emph{t} affects recruitment in year \emph{t + 1}, set \code{lag = 1}.}
    \item{$process_model}{Process model for the ecov time-series. \code{"rw"} = random walk, \code{"ar1"} = 1st order autoregressive, \code{NA} = do not fit}
    \item{$where}{Where does the ecov affect the population? \code{"recuit"} = recruitment,
    \code{"M"} = natural mortality, \code{"growth"} = growth.}
    \item{$how}{How does the ecov affect the \code{$where} process? These options are
    specific to the \code{$where} process.}
    \item{$link_model}{Model describing ecov effect on the \code{$where} process. Options: 'linear' (default) or 'poly-x'
    where x = 2, ... (e.g. 'poly-2' specifies a quadratic model, \eqn{b0 + b1*ecov + b2*ecov^2 + ...}).}
  }

\code{ecov$how} specifies HOW the ecov affects the \code{ecov$where} process.
  \describe{
    \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
    \item{= 1}{"controlling" (dens-indep mortality)}
    \item{= 2}{"limiting" (carrying capacity, e.g. ecov determines amount of suitable habitat)}
    \item{= 3}{"lethal" (threshold, i.e. R --> 0 at some ecov value)}
    \item{= 4}{"masking" (metabolic/growth, decreases dR/dS)}
    \item{= 5}{"directive" (e.g. behavioral)}
  }
Options for M:
  \describe{
    \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
    \item{= 1}{effect on mean M (shared across ages)}
  }

\code{selectivity} specifies options for selectivity, to overwrite default options.
If \code{NULL}, logistic selectivity is assumed. If specifed, \code{selectivity} must be a list with the following entries:
  \describe{
    \item{$model}{Selectivity model for each block. Vector with length = number of selectivity blocks. Each entry must be one of: "age-specific", "logistic", "double-logistic", or "decreasing-logistic".}
    \item{$re}{Time-varying (random effects) for each block. Vector with length = number of selectivity blocks.
                 If \code{NULL}, selectivity parameters in all blocks are constant over time and uncorrelated.
                 Each entry of \code{selectivity$re} must be one of the following options, where the selectivity parameters are:
                 \describe{
                   \item{"none"}{(default) are constant and uncorrelated}
                   \item{"iid"}{vary by year and age/par, but uncorrelated}
                   \item{"ar1"}{correlated by age/par (AR1), but not year}
                   \item{"ar1_y"}{correlated by year (AR1), but not age/par}
                   \item{"2dar1"}{correlated by year and age/par (2D AR1)}
                 }
                }
    \item{$initial_pars}{Initial parameter values for each block. List of length = number of selectivity blocks. Each entry must be a vector of length # parameters in the block, i.e. \code{c(2,0.2)} for logistic or \code{c(0.5,0.5,0.5,1,1,0.5)} for age-specific with 6 ages.}
    \item{$fix_pars}{Which parameters to fix at initial values. List of length = number of selectivity blocks. E.g. model with 3 age-specific blocks and 6 ages, \code{list(c(4,5),4,c(2,3,4))} will fix ages 4 and 5 in block 1, age 4 in block 2, and ages 2, 3, and 4 in block 3.}
  }

\code{M} specifies options for natural mortality and can overwrite default options.
If \code{NULL}, the M is defined as fixed and constanct acros time, with age-specifc values defined in \code{basic_info}. To parameterize M-at-age
shared/mirrored among some but not all ages, modify \code{basic_info$map$M_a} after calling \code{prepare_wham_om_input}
(see vignette on natural mortality for more details). If specified, \code{M} must be a list with the following entries:
  \describe{
    \item{$model}{Natural mortality model options are:
                   \describe{
                     \item{"constant"}{(default) A single mean M shared across all ages}
                     \item{"age-specific"}{M unique for each age}
                     \item{"weight-at-age"}{specifies M as a function of weight-at-age, \eqn{M_y,a = exp(b0 + b1*log(W_y,a))}, as in
                       \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)} and
                       \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035}{Miller & Hyun (2018)}.}
                   }
                 }
    \item{$re}{Time- and age-varying (random effects) on M. Note that random effects can only be estimated using \code{\link{fit_wham}} if
               mean M-at-age parameters are (\code{$est_ages} is not \code{NULL}).
                 \describe{
                   \item{"none"}{(default) M constant in time and across ages.}
                   \item{"iid"}{M varies by year and age, but uncorrelated.}
                   \item{"ar1_a"}{M correlated by age (AR1), constant in time.}
                   \item{"ar1_y"}{M correlated by year (AR1), constant all ages.}
                   \item{"2dar1"}{M correlated by year and age (2D AR1), as in \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047}{Cadigan (2016)}.}
                 }
                }
    \item{$initial_means}{Initial/mean M-at-age vector, with length = number of ages (if \code{$model = "age-specific"})
                         or 1 (if \code{$model = "weight-at-age" or "constant"}). If \code{NULL}, initial mean M-at-age values are taken
                         from \code{basic_info}.}
    \item{$est_ages}{Vector of ages to estimate M when using with \code{\link{fit_wham}} (others will be fixed at initial values). E.g. in a model with 6 ages,
                     \code{$est_ages = 5:6} will estimate M for the 5th and 6th ages, and fix M for ages 1-4. If \code{NULL},
                     M at all ages is fixed at \code{M$initial_means} (if not \code{NULL}) or row 1 of the MAA matrix from the ASAP file (if \code{M$initial_means = NULL}).}
    \item{$logb_prior}{(Only if \code{$model = "weight-at-age"}) TRUE or FALSE (default), should a N(0.305, 0.08) prior be
                       used on log_b? Based on Fig. 1 and Table 1 (marine fish) in \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)}.}
  }

\code{NAA_re} specifies options for random effects on numbers-at-age (NAA, i.e. state-space model or not)
If \code{NULL}, a traditional statistical catch-at-age model is fit (NAA = pred_NAA for all ages, deterministic).
To fit a state-space model, specify \code{NAA_re}, a list with the following entries:
  \describe{
    \item{$sigma}{Which ages allow deviations from pred_NAA? Common options are specified with the strings:
                   \describe{
                     \item{"rec"}{Random effects on recruitment (deviations), all other ages deterministic}
                     \item{"rec+1"}{"Full state space" model with 2 values for \code{sigma_a}, one for recruitment and one shared among other ages}
                   }
                  Alternatively, you can specify a more complex structure by entering a vector with length = n.ages, where each entry points to the 
                  NAA_sigma to use for that age. E.g. c(1,2,2,3,3,3) will parameterize 3 values for \code{sigma_a}, with recruitment (age-1) deviations having their
                  own \code{sigma_R}, ages 2-3 sharing \code{sigma_2}, and ages 4-6 sharing \code{sigma_3}.
                 }
    \item{$cor}{Correlation structure for the NAA deviations. Options are:
                 \describe{
                   \item{"iid"}{NAA deviations vary by year and age, but uncorrelated.}
                   \item{"ar1_a"}{NAA deviations correlated by age (AR1).}
                   \item{"ar1_y"}{NAA deviations correlated by year (AR1).}
                   \item{"2dar1"}{NAA deviations correlated by year and age (2D AR1).}
                 }
               }
  }
}
\examples{
\dontrun{
 input = prepare_wham_om_input()
 mod = fit_wham(input, do.fit = FALSE)
 newdata = mod$simulate(complete=TRUE)
}

}
\seealso{
\code{\link{read_asap3_dat}}, \code{\link{fit_wham}}, \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP}, \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}
}
