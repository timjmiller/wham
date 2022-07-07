make_basic_info <- function(base_years = 1982:2021, ages = 1:10, Fhist = "updown", n_feedback_years = 0) { #changed years
    info <- list()
    info$ages <- ages
    info$years <- as.integer(base_years[1] - 1 + 1:(length(base_years) + n_feedback_years))
    info$n_fleets <- 1 
    info$n_indices <- 1
    na <- length(info$ages)
    ny <- length(info$years)
    
    nby <- length(base_years)
    mid <- floor(nby/2)
    #up then down
    if(Fhist == "updown") info$F <- matrix(0.2 + c(seq(0,0.4,length.out = mid),seq(0.4,0,length.out=nby-mid)),nby, info$n_fleets)

    #down then up
    if(Fhist == "downup") info$F <- matrix(0.2 + c(seq(0.4,0,length.out = mid),seq(0,0.4,length.out=nby-mid)),nby, info$n_fleets)

    if(n_feedback_years>0) info$F <- rbind(info$F, info$F[rep(nby, n_feedback_years),, drop = F]) #same F as terminal year for feedback period

    info$catch_cv <- matrix(0.1, ny, info$n_fleets)
    info$catch_Neff <- matrix(200, ny, info$n_fleets)
    
    info$index_cv <- matrix(0.3, ny, info$n_indices)
    info$index_Neff <- matrix(100, ny, info$n_indices)
    info$fracyr_indices <- matrix(0.5, ny, info$n_indices)
    info$index_units <- rep(1, length(info$n_indices)) #biomass
    info$index_paa_units <- rep(2, length(info$n_indices)) #abundance
    
    info$maturity <- t(matrix(1/(1 + exp(-1*(1:na - na/2))), na, ny))

    L <- 100*(1-exp(-0.3*(1:na - 0)))
    W <- exp(-11)*L^3
    nwaa <- info$n_indices + info$n_fleets + 2
    info$waa <- array(NA, dim = c(nwaa, ny, na))
    for(i in 1:nwaa) info$waa[i,,] <- t(matrix(W, na, ny))

    info$fracyr_SSB <- rep(0.25,ny)
    info$q <- rep(0.3, info$n_indices)

    info$selblock_pointer_fleets <- t(matrix(1:info$n_fleets, info$n_fleets, ny))
    info$selblock_pointer_indices <- t(matrix(info$n_fleets + 1:info$n_indices, info$n_indices, ny))
    return(info)
}


set_NAA_pars <- function(pars, input)
{
    if("NAA_sigma" %in% names(pars)) input$par$log_NAA_sigma[] = log(pars$NAA_sigma)
    if("NAA_rho" %in% names(pars)) input$par$trans_NAA_rho[] = wham:::gen.logit(pars$NAA_rho, -1,1, s = 2)
    if("rec_pars" %in% names(pars)) input$par$mean_rec_pars[] = log(pars$rec_pars)
    return(input)
}

set_M_pars <- function(pars, input)
{
    if("NAA_sigma" %in% names(pars)) input$par$log_NAA_sigma[] = log(pars$NAA_sigma)
    if("rec_pars" %in% names(pars)) input$par$mean_rec_pars[] = log(pars$rec_pars)
    return(input)
}

set_Ecov_beta <- function(beta, input)
{
    ind = which(!is.na(input$map$Ecov_beta))
    input$par$Ecov_beta[ind] = beta
    return(input)
}

get_F_from_catch <- function(om, year, catch, Finit = 0.1, maxF = 10){
    rep = om$report()
    #print(year)
    naa = rep$NAA[year,]
    Maa = rep$MAA[year,]
    sel_tot = rep$FAA_tot[year,]/max(rep$FAA_tot[year,])
    waa = om$input$data$waa[om$input$data$waa_pointer_totcatch, year,]
    get_catch = function(log_F, naa, sel, waa, Maa){
        Faa = exp(log_F) * sel_tot
        Zaa = Maa + Faa
        Catch = 0
        for(a  in 1:length(naa)) Catch = Catch + waa[a] * naa[a] * Faa[a] *(1 - exp(-Zaa[a]))/Zaa[a];
        return(Catch)
    }
    obj = function(log_F) (catch - get_catch(log_F, naa, sel_tot, waa, Maa))^2
    opt = try(nlminb(log(Finit), obj))
    if(!is.character(opt)) Fsolve = exp(opt$par)[1] else Fsolve = maxF
    if(Fsolve>10) Fsolve = maxF
    print(paste0("Fsolve: ", Fsolve))
    return(Fsolve)
}

update_om_F = function(om, year, catch){
    rep = om$report()
    year_ind = which(om$years == year)
    Fsolve = get_F_from_catch(om, year_ind, catch)
    #have to be careful if more than one fleet
    sel = rep$FAA[year_ind,,] #n_fleets x n_ages
    sel_all = sel/max(rep$FAA_tot[year_ind,])
    FAA_all = Fsolve * rbind(sel_all)
    F_fleet = apply(FAA_all, 1, max)
    if(year_ind>1) om$input$par$F_devs[year_ind-1,] = log(F_fleet) - log(rep$F[year_ind-1])
    else om$input$par$log_F1 = log(F_fleet)
    om <- fit_wham(om$input, do.fit = FALSE, MakeADFun.silent = TRUE)

}