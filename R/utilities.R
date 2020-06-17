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
#' under review, 2020.
#'
#' @param X a data matrix
#' @param align_order a bloolean variable, indicating whether to adjust the order of magnitude of parameters.
#'                    NOTE: It is specially designed for \code{\link{MVSKtilting}} function.
#'                    
#' @return A list containing the following elements:
#' \item{\code{mu}}{mean vector.}
#' \item{\code{Sgm}}{Covariance matrix.}
#' \item{\code{Phi}}{Coskewness matrix in vector form (collecting only the unique elements).}
#' \item{\code{Psi}}{Cokurtosis matrix in vector form (collecting only the unique elements).}
#' \item{\code{Phi_shred}}{Partition on \code{Phi}, see reference.}
#' \item{\code{Psi_shred}}{Partition on \code{Psi}, see reference.}
#'
#'
#' @examples
#' \dontrun{
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' params <- estimate_moments_MVSK(X50)
#' }
#'
#' @import PerformanceAnalytics
#' @importFrom stats cov
#' @export
estMomParams <- function(X, align_order = FALSE) {
  N <- ncol(X)
  
  mu  <- colMeans(X)
  Sgm <- cov(X)
  Phi <- PerformanceAnalytics::M3.MM(X, as.mat = FALSE)
  Psi <- PerformanceAnalytics::M4.MM(X, as.mat = FALSE)
  
  if (align_order) {
    d <- abs(evalMoms(w = rep(1/N, N), mom_params = estMomParams(X, FALSE)))
    mu  <- mu  / d[1]
    Sgm <- Sgm / d[2]
    Phi <- Phi / d[3]
    Psi <- Psi / d[4]
  }
  
  Phi_mat <- PerformanceAnalytics:::M3.vec2mat(Phi, N)
  Phi_shred <- lapply(1:N, function(i) Phi_mat[, (1:N)+N*(i-1)])
  
  Psi_mat <- PerformanceAnalytics:::M4.vec2mat(Psi, N)
  Psi_shred <- list()
  for (i in 1:N) {
    tmp <- Psi_mat[, (1:N^2)+N^2*(i-1)]
    Psi_shred[[i]] <- PerformanceAnalytics:::M3.mat2vec(tmp)
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
#' under review, 2020.
#'
#' @param w a numerical vector as portfolio weights.
#' @param mom_params a list of moment parameters, see \code{\link{estMomParams}}.
#'                    
#' @return Four moments of the given portfolio.
#' 
#' @examples
#' \dontrun{
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' params <- estimate_moments_MVSK(X50)
#' moments <- evalMoms(rep(1/50, 50), params)
#' }
#'
#' @import PerformanceAnalytics
#' @export
evalMoms <- function(w, mom_params) {
  c(mean     = as.numeric(w %*% mom_params$mu),
    variance = as.numeric(w %*% mom_params$Sgm %*% w),
    skewness = as.numeric(PerformanceAnalytics:::portm3(w = w, M3 = mom_params$Phi)),
    kurtosis = as.numeric(PerformanceAnalytics:::portm4(w = w, M4 = mom_params$Psi)))
}