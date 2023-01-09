#' Estimate the thresholds from the empirical distributions
#' @param x an ordinal vector
#' @return a vector of thresholds
#' @export

find_thresholds <- function(x) {
  pt <- prop.table(table(x))
  names(pt) <- NULL
  return(c(-20, stats::qnorm(cumsum(pt[-length(pt)])), 20))
}

