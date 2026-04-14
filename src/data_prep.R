#' Rescale Vector to [-0.5, 0.5] Range
#' @param v Numeric vector to rescale
#' @return Rescaled vector with range [-0.5, 0.5]
rescale_to_unit <- function(v) {
  v_min <- min(v, na.rm = TRUE)
  v_max <- max(v, na.rm = TRUE)
  (v - v_min) / (v_max - v_min) - 0.5
}
