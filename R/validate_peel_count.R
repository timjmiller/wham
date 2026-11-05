validate_peel_count <- function(peel, years) {
  if (
    length(peel) != 1L ||
      !is.numeric(peel) ||
      !is.finite(peel) ||
      peel != floor(peel) ||
      peel < 1L ||
      peel > length(years) - 2L
  ) {
    stop("peel must be a positive integer retaining at least two model years.")
  }
}
