# Normalize legacy WHAM output to the current comparison format.
normalize_wham_fit <- function(x) {
  if (length(dim(x$MAA)) == 2L) {
    # Pre-spatial fits represent a single stock and region as year-age matrices.
    x$MAA <- array(x$MAA, dim = c(1L, 1L, dim(x$MAA)))
    for (nm in c("NAA", "log_NAA_rep")) {
      if (!is.null(x[[nm]])) {
        expand <- function(value) array(value, dim = c(1L, 1L, dim(value)))
        x[[nm]]$est <- expand(x[[nm]]$est)
        x[[nm]]$se <- expand(x[[nm]]$se)
        x[[nm]]$ci <- lapply(x[[nm]]$ci, expand)
      }
    }
    if (is.null(x[["log_SSB_all"]])) {
      x$log_SSB_all <- x[["log_SSB"]]
    }
    if (is.null(x[["log_F_tot"]])) {
      x$log_F_tot <- x[["log_F"]]
    }
    for (nm in c("log_SSB_all", "log_F_tot")) {
      if (!is.null(x[[nm]])) {
        x[[nm]]$est <- as.vector(x[[nm]]$est)
        x[[nm]]$se <- as.vector(x[[nm]]$se)
        x[[nm]]$ci <- lapply(x[[nm]]$ci, as.vector)
      }
    }
  }
  # Fits without sdreport may only report SSB and NAA on the natural scale.
  # Normalize these here so every comparison panel sees the same output.
  log_report <- function(est) {
    se <- est
    se[] <- 0
    list(est = log(est), se = se, ci = list(lo = se, hi = se))
  }
  if (is.null(x[["log_SSB_all"]]) && !is.null(x[["SSB"]])) {
    x$log_SSB_all <- log_report(rowSums(cbind(x$SSB$est)))
  }
  if (is.null(x[["log_NAA_rep"]]) && !is.null(x[["NAA"]])) {
    x$log_NAA_rep <- log_report(x$NAA$est)
  }
  x
}

# Find report indices, using a legacy report name when necessary.
wham_report_index <- function(summary, current, legacy) {
  index <- which(rownames(summary) == current)
  if (!length(index)) {
    index <- which(rownames(summary) == legacy)
  }
  index
}
