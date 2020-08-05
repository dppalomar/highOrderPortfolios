# This file is the collection of utilities used in this package

#' @title Estimate first four moment parameters of multivariate observations
#'
#' @description Estimate first four moments of multivariate observations, namely,
#' mean vector, covariance matrix, coskewness matrix, and cokurtosis matrix.
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' Rui Zhou and and Daniel P. Palomar, “Solving High-Order Portfolios via Successive Convex Approximation Algorithms,”
#' under review, 2020. <https://arxiv.org/abs/2008.00863>
#'
#' @param X Data matrix.
#' @param adjust_magnitude Boolean indicating whether to adjust the order of magnitude of parameters.
#'                         NOTE: It is specially designed for \code{\link{design_MVSKtilting_portfolio}} function.
#'                    
#' @return A list containing the following elements:
#' \item{\code{mu}}{Mean vector.}
#' \item{\code{Sgm}}{Covariance matrix.}
#' \item{\code{Phi}}{Co-skewness matrix in vector form (collecting only the unique elements).}
#' \item{\code{Psi}}{Co-kurtosis matrix in vector form (collecting only the unique elements).}
#' \item{\code{Phi_shred}}{Partition on \code{Phi} (see reference).}
#' \item{\code{Psi_shred}}{Partition on \code{Psi} (see reference).}
#'
#'
#' @examples
#' \donttest{
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' X_moments <- estimate_moments(X50)
#' }
#'
#' @importFrom PerformanceAnalytics M3.MM
#' @importFrom PerformanceAnalytics M4.MM
#' @importFrom stats cov
#' @export
estimate_moments <- function(X, adjust_magnitude = FALSE) {
  N <- ncol(X)
  
  mu  <- colMeans(X)
  Sgm <- cov(X)
  Phi <- M3.MM(X, as.mat = FALSE)
  Psi <- M4.MM(X, as.mat = FALSE)
  
  if (adjust_magnitude) {
    d <- abs(eval_portfolio_moments(w = rep(1/N, N), X_moments = estimate_moments(X, FALSE)))
    mu  <- mu  / d[1]
    Sgm <- Sgm / d[2]
    Phi <- Phi / d[3]
    Psi <- Psi / d[4]
  }
  
  Phi_mat <- PerformanceAnalytics_M3.vec2mat(Phi, N)
  Phi_shred <- lapply(1:N, function(i) Phi_mat[, (1:N)+N*(i-1)])
  
  Psi_mat <- PerformanceAnalytics_M4.vec2mat(Psi, N)
  Psi_shred <- list()
  for (i in 1:N) {
    tmp <- Psi_mat[, (1:N^2)+N^2*(i-1)]
    Psi_shred[[i]] <- PerformanceAnalytics_M3.mat2vec(tmp)
    gc()
  }
  gc()
  
  return(list(mu = mu, Sgm = Sgm, Phi = Phi, Psi = Psi, Phi_shred = Phi_shred, Psi_shred = Psi_shred))
}




#' @title Evaluate first four moments of a given portfolio
#'
#' @description Estimate first four moments of a given portfolio's return, namely,
#' mean, variance, skewness, and kurtosis.
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' Rui Zhou and and Daniel P. Palomar, “Solving High-Order Portfolios via Successive Convex Approximation Algorithms,”
#' under review, 2020. <https://arxiv.org/abs/2008.00863>
#'
#' @param w Numerical vector with portfolio weights.
#' @param X_moments List of moment parameters, as obtained with function \code{\link{estimate_moments}}.
#'                    
#' @return Four moments of the given portfolio.
#' 
#' @examples
#' \donttest{
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' X_moments <- estimate_moments(X50)
#' w_moments <- eval_portfolio_moments(rep(1/50, 50), X_moments)
#' }
#'
#' @export
eval_portfolio_moments <- function(w, X_moments) {
  c(mean     = as.numeric(w %*% X_moments$mu),
    variance = as.numeric(w %*% X_moments$Sgm %*% w),
    skewness = as.numeric(PerformanceAnalytics_portm3(w = w, M3 = X_moments$Phi)),
    kurtosis = as.numeric(PerformanceAnalytics_portm4(w = w, M4 = X_moments$Psi)))
}