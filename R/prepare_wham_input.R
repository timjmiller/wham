#' Prepare input data and parameters for WHAM model
#'
#' After the data file is read into R by \code{\link{read_asap3_dat}}, this function
#' prepares the data and parameter settings for \code{\link{fit_wham}}.
#' By default, this will set up a SCAA version like \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP}.
#' As of version 1.0.5, if an asap3 object is not provided, a dummy input is generated with some arbitrary
#' assumptions. The various options described below, such as \code{NAA_re} and \code{selectivity},
#' can still be used without the asap3 object.
#'
#' \code{recruit_model} specifies the stock-recruit model. See \code{wham.cpp} for implementation.
#'   \describe{
#'     \item{= 1}{SCAA (without NAA_re option specified) or Random walk (if NAA_re$sigma specified), i.e. predicted recruitment in year i = recruitment in year i-1}
#'     \item{= 2}{(default) Random about mean, i.e. steepness = 1}
#'     \item{= 3}{Beverton-Holt}
#'     \item{= 4}{Ricker}
#'   }
#' Note: we allow fitting a SCAA model (\code{NAA_re = NULL}), which estimates recruitment in every year as separate fixed effect parameters,
#' but in that case no stock-recruit function is estimated. A warning message is printed if \code{recruit_model > 2} and \code{NAA_re = NULL}.
#' If you wish to use a stock-recruit function for expected recruitment, choose recruitment deviations as random effects,
#' either only age-1 (\code{NAA_re = list(sigma="rec")}) or all ages (\code{NAA_re = list(sigma="rec+1")}, "full state-space" model).
#' See below for details on \code{NAA_re} specification.
#'
#' \code{ecov} specifies any environmental covariate data and model. Environmental covariate data need not span
#' the same years as the fisheries data. It can be \code{NULL} if no environmental data are to be fit.
#' Otherwise, it must be a named list with the following components:
#'   \describe{
#'     \item{$label}{Name(s) of the environmental covariate(s). Used in printing.}
#'     \item{$mean}{Mean observations (matrix). number of years x number of covariates. Missing values = NA.}
#'     \item{$logsigma}{Configure observation standard errors. Options:
#'       \describe{
#'         \item{Matrix of \eqn{log} standard errors with same dimensions as \code{$mean}}{Specified values for each time step }
#'         \item{log standard errors for each covariate, numeric vector or matrix w/ dim 1 x n.ecov}{Specified value the same for all time steps}
#'         \item{estimation option (for all covariates). character string:}{
#'           \code{"est_1"}: Estimated, one value shared among time steps.
#'           \code{"est_re"}: Estimated value for each time step as random effects with two parameters (mean, var)}
#'         \item{list of two elements.}{
#'           First is the matrix of log standard errors or the vector of single values for each covariate as above. 
#'           Second is a character vector of estimation options (\code{NA}, \code{"est_1"},\code{"est_re"}) for each covariate. 
#'           For covariates with non-NA values, values in the first element are ignored.}
#'       }
#'     }
#'     \item{$year}{Years corresponding to observations (vector of same length as \code{$mean} and \code{$logsigma})}
#'     \item{$use_obs}{T/F (or 0/1) vector/matrix of the same dimension as \code{$mean} and \code{$logsigma}.
#'     Use the observation? Can be used to ignore subsets of the ecov without changing data files.}
#'     \item{$lag}{integer vector of offsets between the ecov observation/process and their affect on the stock.
#'     I.e. if SST in year \emph{t} affects recruitment in year \emph{t + 1}, set \code{lag = 1}. May also be a list (length=n_Ecov) of vectors (length = 2+n_indices) if multiple effects of one or more Ecovs are modeled.}
#'     \item{$process_model}{Process model for the ecov time-series. \code{"rw"} = random walk, \code{"ar1"} = 1st order autoregressive, \code{NA} = do not fit}
#'     \item{$where}{Character vector for where each ecov affects the population. \code{"recruit"} = recruitment,
#'     \code{"M"} = natural mortality, \code{"q"} = catchability for indices, \code{"none"} = fit ecov process model(s) but without an effect on the
#'     population. May also be a list (element for each ecov) of character vectors ("none", "recruit", "M", and/or "q") so each ecov can multiple effects.}
#'     \item{$indices}{indices that each ecov affects. Must be a list (length = n_Ecov), where each element is a vector of indices (1:n_indices). Must be provided when any of \code{where} = "q"}
#'     \item{$link_model}{vector of (orthogonal polynomial order) models for effects of each ecov on the \code{$where} process. Options: "none", "linear" (default) or "poly-x"
#'     where x = 2, ... (e.g. "poly-2" specifies a quadratic model, \eqn{b_0 + b_1*ecov + b_2*ecov^2 + ...}). Or a list (length = n_Ecov) of character vectors (same options) for modeling
#'      effects on (1) recruitment, (2) M, and catchabilities for (3) index 1,..., (2+n_indices) index n_indices (length of the vector is 2 + n_indices).}
#'     \item{$ages}{Ages that each ecov affects. Must be a list of length n.ecov, where each element is a vector of ages. Default = list(1:n.ages) to affect all ages.}
#'     \item{$how}{How does the ecov affect the \code{$where} process? An integer vector (length = n_Ecov). If corresponding \code{$where} = "none", then this is ignored.
#'     These options are specific to the \code{$where} process.
#'       \describe{
#'         \item{Recruitment options (see \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}):}{
#'           \describe{
#'             \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
#'             \item{= 1}{"controlling" (dens-indep mortality)}
#'             \item{= 2}{"limiting" (carrying capacity, e.g. ecov determines amount of suitable habitat)}
#'             \item{= 3}{"lethal" (threshold, i.e. R --> 0 at some ecov value)}
#'             \item{= 4}{"masking" (metabolic/growth, decreases dR/dS)}
#'             \item{= 5}{"directive" (e.g. behavioral)}
#'           }}
#'         \item{M options:}{
#'           \describe{
#'             \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
#'             \item{= 1}{effect on mean M (shared across ages)}
#'           }}
#'         \item{Catchability options:}{
#'           \describe{
#'             \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
#'             \item{= 1}{effect on one or more catchabilities (see \code{$indices)})}
#'           }}
#'       }
#'     }
#'   }
#'
#' \code{selectivity} specifies options for selectivity, to overwrite existing options specified in the ASAP data file.
#' If \code{NULL}, selectivity options from the ASAP data file are used. \code{selectivity} is a list with the following entries:
#'   \describe{
#'     \item{$model}{Selectivity model for each block. Vector with length = number of selectivity blocks. Each entry must be one of: "age-specific", "logistic", "double-logistic", or "decreasing-logistic".}
#'     \item{$re}{Time-varying (random effects) for each block. Vector with length = number of selectivity blocks.
#'                  If \code{NULL}, selectivity parameters in all blocks are constant over time and uncorrelated.
#'                  Each entry of \code{selectivity$re} must be one of the following options, where the selectivity parameters are:
#'                  \describe{
#'                    \item{"none"}{(default) are constant and uncorrelated}
#'                    \item{"iid"}{vary by year and age/par, but uncorrelated}
#'                    \item{"ar1"}{correlated by age/par (AR1), but not year}
#'                    \item{"ar1_y"}{correlated by year (AR1), but not age/par}
#'                    \item{"2dar1"}{correlated by year and age/par (2D AR1)}
#'                  }
#'                 }
#'     \item{$initial_pars}{Initial parameter values for each block. List of length = number of selectivity blocks. Each entry must be a vector of length # parameters in the block, i.e. \code{c(2,0.2)} for logistic or \code{c(0.5,0.5,0.5,1,1,0.5)} for age-specific with 6 ages. Default is to set at middle of parameter range. This is 0.5 for age-specific and n.ages/2 for logistic, double-logistic, and decreasing-logistic.}
#'     \item{$fix_pars}{Which parameters to fix at initial values. List of length = number of selectivity blocks. E.g. model with 3 age-specific blocks and 6 ages, \code{list(c(4,5),4,c(2,3,4))} will fix ages 4 and 5 in block 1, age 4 in block 2, and ages 2, 3, and 4 in block 3.}
#'     \item{$n_selblocks}{How many selectivity blocks. Optional. If unspecified and no asap3 object, then this is set to the number of fleets + indices. If specified, ensure other components of \code{selectivity} are consistent.}
#'   }
#'
#' \code{M} specifies estimation options for natural mortality and can overwrite M-at-age values specified in the ASAP data file.
#' If \code{NULL}, the M-at-age matrix from the ASAP data file is used (M fixed, not estimated). To estimate M-at-age
#' shared/mirrored among some but not all ages, modify \code{input$map$M_a} after calling \code{prepare_wham_input}
#' (see vignette for more details). \code{M} is a list with the following entries:
#'   \describe{
#'     \item{$model}{Natural mortality model options are:
#'                    \describe{
#'                      \item{"constant"}{(default) estimate a single mean M shared across all ages}
#'                      \item{"age-specific"}{estimate M_a independent for each age}
#'                      \item{"weight-at-age"}{specifies M as a function of weight-at-age, \eqn{M_y,a = exp(b0 + b1*log(W_y,a))}, as in
#'                        \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)} and
#'                        \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035}{Miller & Hyun (2018)}.}
#'                    }
#'                  }
#'     \item{$re}{Time- and age-varying (random effects) on M. Note that random effects can only be estimated if
#'                mean M-at-age parameters are (\code{$est_ages} is not \code{NULL}).
#'                 \describe{
#'                   \item{"none"}{(default) M constant in time and across ages.}
#'                   \item{"iid"}{M varies by year and age, but uncorrelated.}
#'                   \item{"ar1_a"}{M correlated by age (AR1), constant in time.}
#'                   \item{"ar1_y"}{M correlated by year (AR1), constant all ages.}
#'                   \item{"2dar1"}{M correlated by year and age (2D AR1), as in \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047}{Cadigan (2016)}.}
#'                 }
#'               }
#'     \item{$initial_means}{Initial/mean M-at-age vector, with length = number of ages (if \code{$model = "age-specific"})
#'                          or 1 (if \code{$model = "weight-at-age" or "constant"}). If \code{NULL}, initial mean M-at-age values are taken
#'                          from the first row of the MAA matrix in the ASAP data file.}
#'     \item{$est_ages}{Vector of ages to estimate M (others will be fixed at initial values). E.g. in a model with 6 ages,
#'                      \code{$est_ages = 5:6} will estimate M for the 5th and 6th ages, and fix M for ages 1-4. If \code{NULL},
#'                      M at all ages is fixed at \code{M$initial_means} (if not \code{NULL}) or row 1 of the MAA matrix from the ASAP file (if \code{M$initial_means = NULL}).}
#'     \item{$logb_prior}{(Only if \code{$model = "weight-at-age"}) TRUE or FALSE (default), should a N(0.305, 0.08) prior be
#'                        used on log_b? Based on Fig. 1 and Table 1 (marine fish) in \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)}.}
#'   }
#' 
#' \code{growth} specifies estimation options for somatic growth.
#' If \code{NULL}, the growth parameters are set by default. \code{growth} is a list with the following entries:
#'   \describe{
#'     \item{$model}{Growth model options are:
#'                    \describe{
#'                      \item{"vB_classic"}{(default) estimate five parameters: growth rate (K), asymptotic length (Linf), length-at-age 1 (L1), 
#' 											coefficient of variation of length-at-age 1 (CV1), and coefficient of variation of length-at- the maximum age (CVA).}
#'                      \item{"stock_LAA"}{input mean length-at-age per year for the stock (not available yet, in development).}
#'                    }
#'                  }
#'     \item{$re}{Year- and cohort-varying (random effects) on growth parameters.
#'                 \describe{
#'                   \item{"none"}{(default) Growth parameter constant over years and cohorts.}
#'                   \item{"iid_y"}{Growth parameter varies by year, but uncorrelated.}
#'                   \item{"iid_c"}{Growth parameter varies by cohort, but uncorrelated.}
#'                   \item{"ar1_y"}{Growth parameter correlated by year (AR1).}
#'                   \item{"ar1_c"}{Growth parameter correlated by cohort (AR1).}
#'                 }
#'               }
#'     \item{$init_vals}{Initial growth parameter values, with length = number of parameters. If \code{NULL}, initial values are assumed
#'                       to be K = 0.2, Linf = 100, L1 = 5, CV1 = 0.1, CVA = 0.1.}
#'     \item{$est_pars}{Vector of parameters to be estimated (others will be fixed at initial values). E.g. \code{growth$est_ages = 1:2} will 
#' 						estimate K and Linf, and fix L1, CV1, and CVA. If \code{NULL}, growth parameters are fixed at \code{growth$init_vals}.}
#'   }
#' 
#' \code{LW} specifies estimation options for potential length-weight relationship.
#' If \code{NULL}, the LW parameters are set by default but not used if \code{waa} is provided. \code{LW} is a list with the following entries:
#'   \describe{
#'     \item{$re}{Year- and cohort-varying (random effects) on LW parameters.
#'                 \describe{
#'                   \item{"none"}{(default) LW parameter constant over years and cohorts.}
#'                   \item{"iid_y"}{LW parameter varies by year, but uncorrelated.}
#'                   \item{"iid_c"}{LW parameter varies by cohort, but uncorrelated.}
#'                   \item{"ar1_y"}{LW parameter correlated by year (AR1).}
#'                   \item{"ar1_c"}{LW parameter correlated by cohort (AR1).}
#'                 }
#'               }
#'     \item{$init_vals}{Initial LW parameter values, with length = number of parameters (2). If \code{NULL}, initial values are assumed
#'                       to be a = 5e-06, b = 3.}
#'     \item{$est_pars}{Vector of parameters to be estimated (others will be fixed at initial values). E.g. \code{$est_ages = 1} will 
#' 						estimate a. If \code{NULL}, LW parameters are fixed at \code{LW$init_vals}.}
#'   }
#'
#' \code{NAA_re} specifies options for random effects on numbers-at-age (NAA, i.e. state-space model or not)
#' If \code{NULL}, a traditional statistical catch-at-age model is fit (NAA = pred_NAA for all ages, deterministic).
#' To fit a state-space model, specify \code{NAA_re}, a list with the following entries:
#'   \describe{
#'     \item{$sigma}{Which ages allow deviations from pred_NAA? Common options are specified with the strings:
#'                    \describe{
#'                      \item{"rec"}{Random effects on recruitment (deviations), all other ages deterministic}
#'                      \item{"rec+1"}{"Full state space" model with 2 estimated \code{sigma_a}, one for recruitment and one shared among other ages}
#'                    }
#'                   Alternatively, you can specify a more complex structure by entering a vector with length = n.ages, where each entry points to the
#'                   NAA_sigma to use for that age. E.g. c(1,2,2,3,3,3) will estimate 3 \code{sigma_a}, with recruitment (age-1) deviations having their
#'                   own \code{sigma_R}, ages 2-3 sharing \code{sigma_2}, and ages 4-6 sharing \code{sigma_3}.
#'                  }
#'     \item{$cor}{Correlation structure for the NAA deviations. Options are:
#'                  \describe{
#'                    \item{"iid"}{NAA deviations vary by year and age, but uncorrelated.}
#'                    \item{"ar1_a"}{NAA deviations correlated by age (AR1).}
#'                    \item{"ar1_y"}{NAA deviations correlated by year (AR1).}
#'                    \item{"2dar1"}{NAA deviations correlated by year and age (2D AR1).}
#'                  }
#'                }
#'   }
#' \code{NAA_re} also can be used to configure initial numbers at age and recruitment models. The optional associated components of \code{NAA_re} are:
#'   \describe{
#'		 \item{$N1_model}{Integer determining which way to model the initial numbers at age:
#'			 \describe{
#'					\item{0}{(default) age-specific fixed effects parameters}
#' 					\item{1}{2 fixed effects parameters: an initial recruitment and an instantaneous fishing mortality rate to generate an equilibruim abundance at age.}
#'			 }
#'		 }
#'     \item{$N1_pars}{if N1_model = 0, then these would be the initial values to use for abundance at age in the first year. If N1_model = 1, This would be the
#'				initial numbers in the first age class and the equilibrium fishing mortality rate generating the rest of the numbers at age in the first year.
#'		 }
#'		 \item{$recruit_model}{Integer determining how to model recruitment. Overrides \code{recruit_model} argument to \code{prepare_wham_input}. Must make sure \code{NAA_re$sigma}, \code{NAA_re$cor}
#'				and \code{ecov} are properly specified.
#'			 \describe{
#'				 	 \item{1}{SCAA, estimating all recruitements as fixed effects or a random walk if NAA_re$sigma specified}
#'				 	 \item{2}{estimating a mean recruitment with yearly recruitements as random effects}
#'				 	 \item{3}{Beverton-Holt stock-recruitment with yearly recruitements as random effects}
#'				 	 \item{4}{Ricker stock-recruitment with yearly recruitements as random effects}
#'			 }
#'		 }
#'		 \item{$use_steepness}{T/F determining whether to use a steepness parameterization for a stock-recruit relationship. Only used if recruit_model>2}.
#'		 \item{$recruit_pars}{vector of initial parameters for recruitment model. If use_steepness=F, parameters are "alpha" and "beta"
#'				otherwise they are steepness and R0.
#'		 }
#'   }
#'
#' \code{catchability} specifies options for catchability. If \code{NULL} and \code{asap3} is not NULL, a single catchability parameter for each index is used with initial values specified in ASAP file. If both are NULL, initial catchabilities for all indices = 0.3.
#' Otherwise, it is a list with the following optional entries:
#'   \describe{
#'     \item{$re}{Time-varying (random effects) for each index. Vector with length = number of indices.
#'                  Each entry of \code{catchability$re} must be one of the following options:
#'                  \describe{
#'                    \item{"none"}{(default) are constant}
#'                    \item{"iid"}{vary by year and age/par, but uncorrelated}
#'                    \item{"ar1"}{correlated by year (AR1)}
#'                  }
#'                 }
#'     \item{$initial_q}{Initial catchabilities for each index. vector length = number of indices. Will override values provided in \code{basic_info$q}.
#'        If \code{basic_info$q} and \code{asap3} are not provided, default q values are 0.3.}
#'     \item{$q_lower}{Lower bound for catchabilities for each index. vector length = number of indices. For indices with NULL components default lower values are 0.}
#'     \item{$q_upper}{Upper bound for catchabilities for each index. vector length = number of indices. For indices with NULL components default lower values are 1000.}
#'     \item{$prior_sd}{vector of NA and standard deviations to use for gaussian prior on logit transform of catchability parameter. Length = number of indices.
#'       Indices with non-NA values will have mean logit q as a random effect with mean determined by logit transform of \code{catchability$initial_q} and
#'       sigma as standard error.}
#'   }
#'
#' \code{age_comp} specifies the age composition models for fleet(s) and indices.
#' If \code{NULL}, the multinomial is used because this was the only option in ASAP.
#' The age composition models available are:
#'   \describe{
#'     \item{\code{"multinomial"}}{Multinomial. This is the default because it was the only option in ASAP. 0 parameters.}
#'     \item{\code{"dir-mult"}}{Dirichlet-multinomial. 1 parameter. Effective sample size is estimated by the model (\href{https://www.ccamlr.org/es/system/files/science_journal_papers/07candy.pdf}{Candy 2008})}
#'     \item{\code{"dirichlet-pool0"}}{Dirichlet, pooling zero observations with adjacent age classes. 1. parameter. See \href{https://www.sciencedirect.com/science/article/abs/pii/S0165783613003093}{Francis 2014} and \href{https://cdnsciencepub.com/doi/abs/10.1139/cjfas-2015-0532}{Albertsen et al. 2016}}
#'     \item{\code{"dirichlet-miss0"}}{}{Dirichlet, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-miss0"}}{Logistic normal, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-ar1-miss0"}}{Logistic normal, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-pool0"}}{Logistic normal, pooling zero observations with adjacent age classes. 1 parameter. See \href{https://doi.org/10.1093/icesjms/fsl024}{Schnute and Haigh (2007)} and \href{https://doi.org/10.1016/j.fishres.2013.12.015}{Francis (2014)}}.
#'     \item{\code{"logistic-normal-01-infl"}}{Zero-or-one inflated logistic normal. Inspired by zero-one inflated beta in \href{https://www.sciencedirect.com/science/article/abs/pii/S0167947311003628}{Ospina and Ferrari (2012)}. 3 parameters. . No OSA residuals.}
#'     \item{\code{"logistic-normal-01-infl-2par"}}{Zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters. No OSA residuals.}
#'   }
#' One-step-ahead residuals will be calculated for all but the last two options when \code{do.osa=TRUE} (Nielsen et al. in prep.). An age composition model needs
#' to be specified for each fleet and index. If you would like all fleets and indices to use the same age composition likelihood, you 
#' can simply specify one of the strings above, i.e. \code{age_comp = "logistic-normal-miss0"}. If you do not want the same
#' age composition model to be used for all fleets and indices, you must specify a named list with the following entries:
#'   \describe{
#'     \item{$fleets}{A vector of the above strings with length = the number of fleets.}
#'     \item{$indices}{A vector of the above strings with length = the number of indices.}
#'   }
#' 
#' \code{len_comp} specifies the length composition models for fleet(s) and indices.
#' If \code{NULL}, the multinomial is used.
#' The length composition models available are:
#'   \describe{
#'     \item{\code{"multinomial"}}{Multinomial. This is the default. 0 parameters.}
#'   }
#' 
#' \code{basic_info} is an optional list of information that can be used to set up the population and types of observations when there is no asap3 file given. Particularly useful for setting
#' up an operating model to simulate population processes and observations. Also can be useful for setting up the structure of assessment model when asap3 has not been used.
#' The current options are:
#'   \describe{
#'     \item{$ages}{integer vector of ages (years) with the last being a plus group}
#' 	   \item{$lengths}{vector of fish total length. \code{n_lengths} will be calculated internally.}
#'     \item{$years}{integer vector of years that the population model spans.}
#'     \item{$n_fleets}{number of fleets.}
#'     \item{$agg_catch}{matrix (length(years) x n_fleets) of annual aggregate catches (biomass) for each fleet.}
#'     \item{$catch_paa}{array (n_fleets x length(years) x n_ages) of each fleet's age composition data (numbers).}
#'     \item{$catch_cv}{matrix (length(years) x n_fleets) of annual CVs for each fleet's aggregate catch observations.}
#'     \item{$catch_Neff}{matrix (length(years) x n_fleets) of annual effective sample sizes for each fleet's age composition observation.}
#'     \item{$use_catch_paa}{0/1 matrix (length(years) x n_fleets) indicated whether to use each fleet's age composition observation.}
#'     \item{$catch_pal}{array (n_fleets x length(years) x n_lengths) of each fleet's length composition data (numbers).}
#'     \item{$catch_NeffL}{matrix (length(years) x n_fleets) of annual effective sample sizes for each fleet's length composition observation.}
#'     \item{$use_catch_pal}{0/1 matrix (length(years) x n_fleets) indicated whether to use each fleet's length composition observation.}
#'     \item{$use_catch_alk}{0/1 array (n_fleets x length(years) x n_lengths) indicated whether to use each fleet's conditional age at length composition.}
#'     \item{$catch_alk_Neff}{array (length(years) x n_fleets x n_lengths) of annual effective sample sizes for each fleet's conditional age at length composition.}
#'     \item{$catch_alk}{array (n_fleets x length(years) x n_lengths x n_ages) for each fleet's conditional age at length composition.}
#'     \item{$selblock_pointer_fleets}{integer matrix (length(years) x n_fleets) indicated which selectivity model to use for each fleet each year. Must be consistent with options to \code{selectivity} option.}
#'     \item{$F}{matrix (length(years) x n_fleets) of annual fishing mortality rates for each fleet to initialize the model.}
#'     \item{$n_indices}{number of indices.}
#'     \item{$agg_indices}{matrix (length(years) x n_indices) of annual aggregate catches (biomass or number) for each fleet.}
#'     \item{$index_paa}{array (n_indices x length(years) x n_ages) of each index's age composition data (biomass or number).}
#'     \item{$index_cv}{matrix (length(years) x n_indices) of annual CVs for each index's aggregate observations.}
#'     \item{$index_Neff}{matrix (length(years) x n_indices) of annual effective sample sizes for each index's age composition observation.}
#'     \item{$units_indices}{1/2 matrix (length =  n_indices) indicated whether indices are in biomass or numbers, respectively.}
#'     \item{$units_index_paa}{1/2 matrix (length = n_indices) indicated whether to use each index's age composition observation are in numbers or biomass.}
#'     \item{$use_index_paa}{0/1 matrix (length(years) x n_indices) indicated whether to use each index's age composition observation.}
#'     \item{$index_pal}{array (n_indices x length(years) x n_lengths) of each index's length composition data (number).}
#'     \item{$index_NeffL}{matrix (length(years) x n_indices) of annual effective sample sizes for each index's length composition observation.}
#'     \item{$units_index_pal}{1/2 matrix (length = n_indices) indicated whether to use each index's length composition observation are in numbers or biomass. Should always be 2.}
#'     \item{$use_index_pal}{0/1 matrix (length(years) x n_indices) indicated whether to use each index's length composition observation.}
#'     \item{$use_index_caal}{0/1 array (n_indices x length(years) x n_lengths) indicated whether to use each index's conditional age at length composition.}
#'     \item{$index_alk_Neff}{array (length(years) x n_indices x n_lengths) of annual effective sample sizes for each index's conditional age at length composition.}
#'     \item{$index_alk}{array (n_indices x length(years) x n_lengths x n_ages) for each index's conditional age at length composition.}
#'     \item{$selblock_pointer_indices}{integer matrix (length(years) x n_indices) indicated which selectivity model to use for each index each year. Must be consistent with options to \code{selectivity} option.}
#'     \item{$fracyr_indices}{matrix (length(years) x n_indices) of annual proportions of the year elapsed when each index is observing the population.}
#'     \item{$waa}{array ((n_fleets + n_indices + 2) x length(years) x length(ages)) of annual weight at at age for each fleet, each index, total catch, and spawning biomass. If not provided, the model will use the \code{LW} parameters.}
#'     \item{$maturity}{matrix (length(years) x length(ages)) of annual maturity at age for estimating spawning biomass.}
#'     \item{$fracyr_SSB}{vector (1 or length(years)) of yearly proportions (0-1) of the year at which to calculate spawning biomass.}
#'     \item{$Fbar_ages}{integer vector of ages to use to average total F at age defining fully selected F for reference points. May not be clearly known until a model is fitted.}
#'     \item{$q}{vector (length(n_indices)) of catchabilities for each of the indices to initialize the model.}
#'     \item{$percentSPR}{(0-100) percentage of unfished spawning biomass per recruit for determining equilibrium fishing mortality reference point}
#'     \item{$percentFXSPR}{(0-100) percentage of SPR-based F to use in projections.}
#'     \item{$percentFMSY}{(0-100) percentage of Fmsy to use in projections.}
#'     \item{$XSPR_R_avg_yrs}{which years to average recruitments for calculating SPR-based SSB reference points. Default is 1:length(years)}
#'     \item{$XSPR_R_opt}{1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4): use average R estimates(predictions).}
#'     \item{$simulate_process_error}{T/F vector (length = 7). When simulating from the model, whether to simulate any process errors for NAA, M, selectivity, Ecov, q, growth, and LW. Only used if applicable.}
#'     \item{$simulate_observation_error}{T/F vector (length = 3). When simulating from the model, whether to simulate  catch, index, and ecov observations.}
#'     \item{$simulate_period}{T/F vector (length = 2). When simulating from the model, whether to simulate base period (model years) and projection period.}
#'     \item{$bias_correct_process}{T/F. Perform bias correction of log-normal random effects for NAA.}
#'     \item{$bias_correct_observation}{T/F. Perform bias correction of log-normal observations.}
#'   }
#' If other arguments to \code{prepare_wham_input} are provided such as \code{selectivity}, \code{M}, \code{age_comp}, and \code{len_comp}, the information provided there
#' must be consistent with \code{basic_info}. For example the dimensions for number of years, ages, fleets, and indices.
#'
#' @param asap3 (optional) list containing data and parameters (output from \code{\link{read_asap3_dat}})
#' @param recruit_model numeric, option to specify stock-recruit model (see details)
#' @param model_name character, name of stock/model
#' @param ecov (optional) named list of environmental covariate data and parameters (see details)
#' @param selectivity (optional) list specifying selectivity options by block: models, initial values, parameters to fix, and random effects (see details)
#' @param M (optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see details)
#' @param growth (optional) list specifying growth options: model, random effects, initial values, and parameters to fix (see details)
#' @param LW (optional) list specifying length-weight relationship options: random effects, initial values, and parameters to fix (see details)
#' @param NAA_re (optional) list specifying options for random effect on numbers-at-age, initial numbers at age, and recruitment model (see details)
#' @param catchability (optional) list specifying options for priors and random effects on catchability (see details)
#' @param age_comp (optional) character or named list, specifies age composition model for fleet(s) and indices (see details)
#' @param len_comp (optional) character or named list, specifies length composition model for fleet(s) and indices (see details)
#' @param basic_info (optional) list of basic population information for use when asap3 is not provided (see details)
#'
#' @return a named list with the following components:
#'   \describe{
#'     \item{\code{data}}{Named list of data, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{par}}{Named list of parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{map}}{Named list defining how to optionally collect and fix parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{random}}{Character vector of parameters to treat as random effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{years}}{Numeric vector of years to fit WHAM model (specified in ASAP3 .dat file)}
#'     \item{\code{ages.lab}}{Character vector of age labels, ending with plus-group (specified in ASAP3 .dat file)}
#'     \item{\code{model_name}}{Character, name of stock/model (specified in call to \code{prepare_wham_input})}
#'   }
#'
#' @seealso \code{\link{read_asap3_dat}}, \code{\link{fit_wham}}, \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP}, \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}
#'
#' @examples
#' \dontrun{
#' asap3 = read_asap3_dat("ex1_SNEMAYT.dat")
#' input = prepare_wham_input(asap3)
#' mod = fit_wham(input)
#'
#' # no ASAP3 file, default parameter values and configuration
#' input = prepare_wham_input()
#' mod = fit_wham(input, fit = FALSE)
#' set.seed(8675309)
#' simdata = mod$simulate(complete=TRUE)
#' input$data = simdata
#' fit = fit_wham(input, do.osa = FALSE)
#' }
#'
#' @export
prepare_wham_input <- function(asap3 = NULL, model_name="WHAM for unnamed stock", recruit_model=2, ecov=NULL, selectivity=NULL, 
	growth=NULL, LW = NULL, M=NULL, NAA_re=NULL, LAA = NULL, catchability=NULL, age_comp=NULL, len_comp = NULL, basic_info = NULL){

	data = list()
	par = list()
	map = list()
	random = character()
	input = list(
	  	data = data,
	  	par = par,
	  	map = map,
	  	random = random,
	  	years = NULL, years_full = NULL, ages.lab = NULL, model_name = model_name, asap3 = asap3)


	#if(!is.null(asap3) & !is.null(growth)) stop('Growth feature does not work with ASAP3 input. Please use the basic_info argument instead.')

	if(is.null(basic_info)) basic_info = list(recruit_model = recruit_model)
	else basic_info$recruit_model = recruit_model

	waa_opts = NULL
	waa_names = c("waa")
	if(any(names(basic_info) %in% waa_names)) waa_opts = basic_info[waa_names]

	catch_opts = NULL
	catch_names = c("n_fleets","agg_catch", "catch_paa", "catch_cv","catch_Neff", "use_catch_paa",
					"catch_pal", "catch_NeffL", "use_catch_pal", "catch_caal", "catch_caal_Neff", "use_catch_caal",
					"selblock_pointer_fleets")
	if(any(names(basic_info) %in% catch_names)) catch_opts = basic_info[catch_names]

	index_opts = NULL
	index_names = c("n_indices", "agg_indices", "index_paa", "fracyr_indices", "index_cv", "index_Neff", "units_indices",
		"units_index_paa", "use_indices", "use_index_paa",
		"index_pal", "index_NeffL", "units_index_pal", "use_index_pal", "index_caal", "index_caal_Neff", "use_index_caal",
		"selblock_pointer_indices")
	if(any(names(basic_info) %in% index_names)) index_opts = basic_info[index_names]

	F_opts = NULL
	F_names = c("F")
	if(any(names(basic_info) %in% F_names)) F_opts = basic_info[F_names]

	q_opts = catchability
	if(any(names(basic_info) == "q") & !any(names(q_opts) == "initial_q")) q_opts$initial_q = basic_info$q

	if(!is.null(asap3))
	{
	  asap3 = asap3$dat
  	  input$asap3 = asap3
	  input$data$n_ages = asap3$n_ages
	  input$data$lengths = seq(from = 2, to = 30, by = 2)
	  input$data$n_lengths = length(input$data$lengths)
	  input$data$fracyr_SSB = rep(asap3$fracyr_spawn, asap3$n_years)
	  input$data$mature = asap3$maturity
	  input$data$Fbar_ages = seq(asap3$Frep_ages[1], asap3$Frep_ages[2])
  	input$years <- asap3$year1 + 1:asap3$n_years - 1
	}
	else
	{
		#if no asap3 is provided, make some default values to
		input = add_basic_info(input, basic_info)
	}

	# print("start")
	#some basic input elements see the function code below
	input = initial_input_fn(input, basic_info)

	# Catch
	input = set_catch(input, catch_opts)
	#print("catch")

	# Indices/surveys
	input = set_indices(input, index_opts)
	#print("indices")

	# WAA in case we want to modify how weight-at age is handled
	input = set_WAA(input, waa_opts)
	#print("WAA")

	# NAA and recruitment options
	input = set_NAA(input, NAA_re)
	#print("NAA")

	input = set_q(input, q_opts)
	#print("q")

	# Selectivity
	input = set_selectivity(input, selectivity)
	#print("selectivity")

	# Growth
	input = set_growth(input, growth, LAA)
	#print("growth")

	# Growth
	input = set_LAA(input, LAA)
	#print("growth")

	# LW
	input = set_LW(input, LW)
	#print("LW")

	# Age composition model
	input = set_age_comp(input, age_comp)
	#print("age_comp")

	# Length composition model
	input = set_len_comp(input, len_comp)
	#print("len_comp")

	#in case we want to add alternative F options
	input = set_F(input, F_opts)
	#print("F")

	#set up natural mortality
	input = set_M(input, M)
	#print("M")

	#set up ecov data and parameters. Probably want to make sure to do this after set_NAA.
	input = set_ecov(input, ecov)
	#print("ecov")

	# add vector of all observations for one step ahead residuals ==========================
	input = set_osa_obs(input)
	#print("osa_obs")

	# projection data will always be modified by 'prepare_projection'
	input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
	#print("proj")

	#set any parameters as random effects
	input = set_random(input)
	#print("random")

	return(input)
}


gen.logit <- function(x, low, upp, s=1) (log((x-low)/(upp-x)))/s


initial_input_fn = function(input, basic_info){
	#this function is a helper so that this code can be run in other set up functions like change_wham_input
  input$years_full = input$years
  input$ages.lab = paste0(1:input$data$n_ages, c(rep("",input$data$n_ages-1),"+"))
	if(!is.null(basic_info$ages)) {
		if(!is.integer(basic_info$ages) | length(basic_info$ages) != input$data$n_ages) stop("basic_info$ages has been specified, but it is not an integer vector or it is not = n_ages")
		else {
  		input$ages.lab = paste0(basic_info$ages, c(rep("",input$data$n_ages-1),"+"))
		}
	}

  input$data$n_years_model = length(input$years)
  input$data$n_years_catch = length(input$years)
  input$data$n_years_indices = length(input$years)
  input$data$recruit_model = basic_info$recruit_model #this is made from argument of the same name to prepare_wham_input
  input$data$which_F_age = rep(input$data$n_ages,input$data$n_years_model) #plus group by default used to define full F (total) IN annual reference points for projections, only. prepare_projection changes it to properly define selectivity for projections.
  input$data$which_F_age_static = input$data$n_ages #plus group by default used to define full F (total) for static SPR-based ref points.

  input$data$bias_correct_pe = 1 #bias correct log-normal process errors?
  input$data$bias_correct_oe = 1 #bias correct log-normal observation errors?
  input$data$simulate_state = rep(1,7) #simulate state variables (NAA, M, sel, Ecov, q, growth (or LAA), LW)
  input$data$simulate_data = rep(1,3) #simulate data types (catch, indices, Ecov)
  input$data$simulate_period = c(1,0) #simulate above items for (model years, projection years)
  input$data$percentSPR = 40 #percentage of unfished SSB/R to use for SPR-based reference points
  input$data$percentFXSPR = 100 # percent of F_XSPR to use for calculating catch in projections
  input$data$percentFMSY = 100 # percent of F_XSPR to use for calculating catch in projections
  # data$XSPR_R_opt = 3 #1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4): use average R estimates(predictions). See next line for years to average over.
  input$data$XSPR_R_opt = 2 # default = use average R estimates
  input$data$XSPR_R_avg_yrs = 1:input$data$n_years_model-1 #model year indices to use for averaging recruitment when defining SSB_XSPR (if XSPR_R_opt = 2,4)
	input$data$static_FXSPR_init = 0.1 #initial value for Newton search of static F (spr-based) reference point (inputs to spr are averages of annual values using avg_years_ind)
  
  if(!is.null(basic_info$bias_correct_process)) input$data$bias_correct_pe = basic_info$bias_correct_process
  if(!is.null(basic_info$bias_correct_observation)) input$data$bias_correct_oe = basic_info$bias_correct_observation
  if(!is.null(basic_info$simulate_process_error)) input$data$simulate_state = basic_info$simulate_process_error
  if(!is.null(basic_info$simulate_observation_error)) input$data$simulate_data = basic_info$simulate_observation_error
  if(!is.null(basic_info$simulate_period)) input$data$simulate_period = basic_info$simulate_period

  if(!is.null(basic_info$percentSPR)) input$data$percentSPR = basic_info$percentSPR
  if(!is.null(basic_info$percentFXSPR)) input$data$percentFXSPR = basic_info$percentFXSPR
  if(!is.null(basic_info$percentFMSY)) input$data$percentFMSY = basic_info$percentFMSY
  if(!is.null(basic_info$XSPR_R_opt)) input$data$XSPR_R_opt = basic_info$XSPR_R_opt
  if(!is.null(basic_info$XSPR_R_avg_yrs)) input$data$XSPR_R_avg_yrs = basic_info$XSPR_R_avg_yrs - 1 #user input shifted to start @ 0 

  return(input)

}
