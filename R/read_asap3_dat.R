#' Read an ASAP3 .dat file into R
#'
#' WHAM is built on ASAP ([Legault and Restrepo (1999)](http://sedarweb.org/docs/wsupp/S12RD06\%20ASAPdoc.pdf)) and this
#' function provides functionality to use a preexisting ASAP3 input data file. The
#' output of \code{read_asap3_dat} should then be passed to \code{\link{prepare_wham_input}}.
#' If you are not familiar with ASAP3 input files, see the ASAP \href{https://github.com/cmlegault/ASAPplots/tree/master/pdf}{documentation}
#' and \href{https://nmfs-fish-tools.github.io/ASAP/}{code}.
#'
#' @param filename character vector, names of ASAP3 .dat files. The file either needs to be
#' in the current working directory, or \code{filename} can include the path. If multipile files,
#' a multi-stock model will be assumed.
#'
#' @return a named list with the following components:
#'   \describe{
#'     \item{\code{dat}}{Named list of input data and parameters}
#'     \item{\code{comments}}{Comments at top of ASAP3 .dat file (indicated by "#")}
#'   }
#'
#' @seealso \code{\link{prepare_wham_input}}, \code{\link{fit_wham}}, \href{https://github.com/cmlegault/ASAPplots/tree/master/pdf}{ASAP documentation}
#'
#' @examples
#' \dontrun{
#' asap3 = read_asap3_dat("ASAP_SNEMAYT.dat")
#' input = prepare_wham_input(asap3)
#' mod = fit_wham(input)
#' }
#'
#' @export
read_asap3_dat <- function(filename){

  single_stock_fun = function(filename){
    char.lines <- readLines(filename)
    com.ind <- which(substring(char.lines,1,1) == "#")
    #print(com.ind)
    dat.start <- com.ind[c(which(diff(com.ind)>1), length(com.ind))]
    comments <- char.lines[dat.start]
    #print(comments)
    #print(dat.start)
    #print(length(dat.start))
    dat <- list()
    ind <- 0
    dat$n_years <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$year1 <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$n_ages <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$n_fleets <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    #print(dat)
    #print(ind)
    dat$n_fleet_sel_blocks <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$n_indices <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)

    dat$M <- matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_years*dat$n_ages), dat$n_years, dat$n_ages, byrow = TRUE)
    dat$fec_opt <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$fracyr_spawn <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$maturity <- matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_years*dat$n_ages), dat$n_years, dat$n_ages, byrow = TRUE)
    dat$n_WAA_mats <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$WAA_mats <- lapply(1:dat$n_WAA_mats, function(x) matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind+x], n = dat$n_years*dat$n_ages), dat$n_years, dat$n_ages, byrow = TRUE))

    ind <- ind+dat$n_WAA_mats
    npt <- dat$n_fleets * 2 + 2 + 2
    dat$WAA_pointers <- sapply(1:npt, function(x) scan(filename, quiet=T, what = integer(), skip = dat.start[ind+1]+x-1, n = 1))
    ind <- ind + 1
    # print(ind)

    dat$sel_block_assign <- lapply(1:dat$n_fleets, function(x) scan(filename, quiet=T, what = integer(), skip = dat.start[ind+x], n = dat$n_years))
    ind <- ind+dat$n_fleets
    dat$sel_block_option <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_fleet_sel_blocks)
    # print(ind)
    # print(dat.start[ind])
    dat$sel_ini <- lapply(1:dat$n_fleet_sel_blocks, function(x) matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind+x], n = 4*(dat$n_ages+6)), dat$n_ages+6, 4, byrow = TRUE))
    ind <- ind + dat$n_fleet_sel_blocks
    dat$fleet_sel_start_age <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)
    dat$fleet_sel_end_age <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)
    dat$Frep_ages <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 2)
    dat$Frep_type <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$use_like_const <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$release_mort <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)

    dat$CAA_mats <- lapply(1:dat$n_fleets, function(x) matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind+x], n = dat$n_years*(dat$n_ages+1)), dat$n_years, dat$n_ages+1, byrow = TRUE))
    ind <- ind + dat$n_fleets
    dat$DAA_mats <- lapply(1:dat$n_fleets, function(x) matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind+x], n = dat$n_years*(dat$n_ages+1)), dat$n_years, dat$n_ages+1, byrow = TRUE))
    ind <- ind + dat$n_fleets
    dat$prop_rel_mats <- lapply(1:dat$n_fleets, function(x) matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind+x], n = dat$n_years*(dat$n_ages)), dat$n_years, dat$n_ages, byrow = TRUE))
    ind <- ind + dat$n_fleets
    # print(ind)

    dat$index_units <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$index_acomp_units <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$index_WAA_pointers <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$index_month <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$index_sel_choice <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$index_sel_option <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$index_sel_start_age <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$index_sel_end_age <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$use_index_acomp <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$use_index <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$index_sel_ini <- lapply(1:dat$n_indices, function(x) matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind+x], n = 4*(dat$n_ages+6)), dat$n_ages+6, 4, byrow = TRUE))
    ind <- ind + dat$n_indices
    # print(dat$n_indices)
    #stop()
    # print(dat$index_sel_ini)
    dat$IAA_mats <- lapply(1:dat$n_indices, function(x) matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind+x], n = dat$n_years*(dat$n_ages+4)), dat$n_years, dat$n_ages+4, byrow = TRUE))
    ind <- ind + dat$n_indices
    # print(ind)

    dat$phase_F1 <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$phase_F_devs <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$phase_rec_devs <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$phase_N1_devs <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$phase_q <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$phase_q_devs <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$phase_SR_scalar <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$phase_steepness <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$recruit_cv <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_years)

    dat$lambda_index <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$lambda_catch <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)
    dat$lambda_discard <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)

    dat$catch_cv <- matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_years*dat$n_fleets), dat$n_years, dat$n_fleets, byrow = TRUE)
    dat$discard_cv <- matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_years*dat$n_fleets), dat$n_years, dat$n_fleets, byrow = TRUE)
    dat$catch_Neff <- matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_years*dat$n_fleets), dat$n_years, dat$n_fleets, byrow = TRUE)
    dat$discard_Neff <- matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_years*dat$n_fleets), dat$n_years, dat$n_fleets, byrow = TRUE)
    # print(ind)

    dat$lambda_F1 <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)
    dat$cv_F1 <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)
    dat$lambda_F_devs <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)
    dat$cv_F_devs <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)

    dat$lambda_N1_devs <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$cv_N1_devs <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$lambda_rec_devs <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)

    dat$lambda_q <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$cv_q <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)

    dat$lambda_q_devs <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$cv_q_devs <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)

    dat$lambda_steepness <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$cv_steepness <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)

    dat$lambda_SR_scalar <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$cv_SR_scalar <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    # print(ind)

    dat$N1_flag <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$N1_ini <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_ages)
    dat$F1_ini <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)
    dat$q_ini <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = dat$n_indices)
    dat$SR_scalar_type <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$SR_scalar_ini <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$steepness_ini <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$Fmax <- scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$ignore_guesses <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    # print(ind)

    dat$do_proj <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$dir_fleet <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$n_fleets)
    dat$nfinalyear <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    n <- dat$nfinalyear-dat$year1-dat$n_years+1
    # print(n)
    # print(ind)
    # print(dat.start[ind])
    if(n>0) dat$proj_ini <- matrix(scan(filename, quiet=T, what = double(), skip = dat.start[ind <- ind + 1], n = n*5), n, 5, byrow = TRUE)
    else dat$proj_ini <- matrix(nrow = 0, ncol = 5)
    dat$doMCMC <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$MCMC_nyear_opt <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$MCMC_nboot <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$MCMC_nthin <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$MCMC_nseed <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$fill_R_opt <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$R_avg_start <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$R_avg_end <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$make_R_file <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    dat$testval <- scan(filename, quiet=T, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
    # print(dat$testval)
    # print(ind)
    fleet.names <- index.names <- NULL
    data.labels <- grep("#$",char.lines, fixed = T, value =T)
    if(length(data.labels) == dat$n_fleets+ dat$n_indices){
      data.labels <- sapply(strsplit(data.labels, "#$", fixed = T), paste, collapse = '')
      dat$fleet.names <- data.labels[1:dat$n_fleets]
      dat$index.names <- data.labels[dat$n_fleets + 1:dat$n_indices]
    }

    return(list(dat = dat, comments = comments, fleet.names = dat$fleet.names, index.names = dat$index.names))
  }
  return(lapply(filename, single_stock_fun))
}
