#' @title Estimate first four moments of multivariate observations
#'
#' @description Estimate first four moments of multivariate observations, namely,
#' mean, variance, skewness, and kurtosis.
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' Rui Zhou and and Daniel P. Palomar, “Solving High-Order Portfolios via Successive Convex Approximation Algorithms,”
#' under review, 2020.
#'
#' @examples
#' library(MVSKportfolio)
#' data(X50)
#' 
#' estimate_moments_MVSK(X50)
#'
#' @import PerformanceAnalytics
#' @export
estimate_moments_MVSK <- function(X) {
  return(list(mu = NULL,
              cov = NULL,
              skewness = NULL,
              kurtosis = NULL))
}



#' @title Design high-order portfolio based on weighted linear combination of first four moments
#' 
#' @description Design high-order portfolio based on weighted linear combination of first four moments
#' (i.e., mean, variance, skewness, and kurtosis):
#' \preformatted{
#'   maximize     w'*mu
#'   subject to   w'*Sigma*w <= var_max, 
#'                ||w||_1 <= leverage,  (w >= 0),
#'                (beta'*w = 0),
#' }
#' where \code{beta} can be a vector or a matrix.
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' Rui Zhou and and Daniel P. Palomar, “Solving High-Order Portfolios via Successive Convex Approximation Algorithms,”
#' under review, 2020.
#'
#' @examples
#' library(MVSKportfolio)
#' data(X50)
#' 
#' moments <- estimate_moments_MVSK(X50)
#' w <- design_MVSK_portfolio(moments)
#'
#' @export
design_MVSK_portfolio <- function() {
  
}


