#' Make one or more selectivity blocks with age-specific parameters
#'
#' @param input list containing data and parameters (output from \code{\link{prepare_wham_input}})
#' @param selblocks numeric, number of age-specific selectivity blocks
#'
#' @return a modified list of data and parameters
#'
#' @examples
#' \dontrun{
#' asap3 = read_asap3_dat("ASAP_SNEMAYT.dat")
#' input = prepare_wham_input(asap3)
#' input = set_age_sel0(input, selblocks = 1:3)
#' mod = fit_wham(input)
#' }
#'
#' @export
set_age_sel0 <- function(input, selblocks){
  temp = input
  temp$map$logit_selpars = as.integer(as.character(temp$map$logit_selpars))
  temp$map$logit_selpars = matrix(temp$map$logit_selpars, temp$data$n_selblocks, temp$data$n_ages + 6)
  temp$map$logit_selpars[selblocks,] = NA
  n.par  = sum(!is.na(temp$map$logit_selpars))
  if(n.par) temp$map$logit_selpars[which(!is.na(temp$map$logit_selpars))] = 1:sum(!is.na(temp$map$logit_selpars))
  temp$data$selblock_models[selblocks] = 1
  for(x in selblocks){
    ind = list(fleets = which(apply(temp$data$selblock_pointer_fleets == x,2,sum) > 0))
    ind$indices = which(apply(temp$data$selblock_pointer_indices == x,2,sum) > 0)
    paa = matrix(nrow = 0, ncol = temp$data$n_ages)
    if(length(ind$fleets)) for(f in ind$fleets)
    {
      y = temp$data$catch_paa[f,which(temp$data$selblock_pointer_fleets[,f] == x & temp$data$use_catch_paa[,f] == 1),]
      paa = rbind(paa,y)
    }
    if(length(ind$indices)) for(i in ind$indices)
    {
      y = temp$data$index_paa[i,which(temp$data$selblock_pointer_indices[,i] == x & temp$data$use_index_paa[,i] == 1),]
      paa = rbind(paa,y)
    }
    #print(dim(paa))
    y = apply(paa,2,sum)
    #print(y)
    temp$par$logit_selpars[x,temp$data$n_ages + 1:6] = Inf
    temp$par$logit_selpars[x,which(y < 1e-5)] = -Inf
    temp$par$logit_selpars[x,which(y >= 1e-5)] = 0
    if(sum(y >= 1e-5)) temp$map$logit_selpars[x,which(y >= 1e-5)] = n.par + 1:sum(y >= 1e-5)
    n.par = n.par + sum(y >= 1e-5)
  }
  # temp$map.n = temp$map$logit_selpars
  temp$map$logit_selpars = factor(temp$map$logit_selpars)
  return(temp)
}
