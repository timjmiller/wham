# Pad observation-sigma REs and maps for added Ecov years; leave other fields alone.
extend_ecov_input <- function(par, map, data, n_years_proj_Ecov) {
  if (any(data$Ecov_obs_sigma_opt == 4L)) {
    sigma_map <- map$Ecov_obs_logsigma_re
    if (is.null(sigma_map)) {
      sigma_map <- factor(seq_along(par$Ecov_obs_logsigma_re))
    }
    sigma_map <- matrix(as.integer(sigma_map), ncol = data$n_Ecov)
    new_sigma_map <- matrix(NA_integer_, n_years_proj_Ecov, data$n_Ecov)
    new_sigma <- matrix(0, n_years_proj_Ecov, data$n_Ecov)
    for (i in which(data$Ecov_obs_sigma_opt == 4L)) {
      new_sigma[, i] <- par$Ecov_obs_sigma_par[1L, i]
      new_sigma_map[, i] <- max(0, sigma_map, new_sigma_map, na.rm = TRUE) +
        seq_len(n_years_proj_Ecov)
    }
    par$Ecov_obs_logsigma_re <- rbind(par$Ecov_obs_logsigma_re, new_sigma)
    map$Ecov_obs_logsigma_re <- factor(rbind(sigma_map, new_sigma_map))
  }
  list(par = par, map = map)
}
