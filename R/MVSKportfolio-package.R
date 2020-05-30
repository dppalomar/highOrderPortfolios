#' MVSKportfolio: Design of High-Order Portfolios via Mean, Variance, Skewness, and Kurtosis
#' 
#' The classical Markowitz's mean-variance portfolio formulation ignores 
#' heavy tails and skewness. High-order portfolios use higher order moments to
#' better characterize the return distribution. Different formulations and fast 
#' algorithms are proposed for high-order portfolios based on the mean, variance, 
#' skewness, and kurtosis.
#'
#'
#' @section Functions:
#' \code{\link{estimate_moments_MVSK}} and \code{\link{design_MVSK_portfolio}}
#'
#' @section Help:
#' For a quick help see the README file:
#' \href{https://github.com/dppalomar/MVSKportfolio/blob/master/README.md}{GitHub-README}.
#'
# For more details see the vignette:
# \href{https://CRAN.R-project.org/package=MVSKportfolio/vignettes/xxx.html}{CRAN-vignette}.
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' Rui Zhou and and Daniel P. Palomar, “Solving High-Order Portfolios via Successive Convex Approximation Algorithms,”
#' under review, 2020.
#'
#' @docType package
#' @name MVSKportfolio-package
NULL
