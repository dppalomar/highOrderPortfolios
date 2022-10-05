#' @title Estimate first four moment parameters of multivariate observations
#'
#' @description Estimate first four moments of multivariate observations, namely,
#' mean vector, covariance matrix, coskewness matrix, and cokurtosis matrix.
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' R. Zhou and D. P. Palomar, "Solving High-Order Portfolios via Successive Convex Approximation Algorithms," 
#' in \emph{IEEE Transactions on Signal Processing}, vol. 69, pp. 892-904, 2021.
#' <doi:10.1109/TSP.2021.3051369>.
#'
#' @param X Data matrix.
#' @param adjust_magnitude Boolean indicating whether to adjust the order of magnitude of parameters.
#'                         Note: this is specially designed for the function \code{\link{design_MVSKtilting_portfolio}()}.
#'                    
#' @return A list containing the following elements:
#' \item{\code{mu}}{Mean vector.}
#' \item{\code{Sgm}}{Covariance matrix.}
#' \item{\code{Phi_mat}}{Co-skewness matrix.}
#' \item{\code{Psi_mat}}{Co-kurtosis matrix.}
#' \item{\code{Phi}}{Co-skewness matrix in vector form (collecting only the unique elements).}
#' \item{\code{Psi}}{Co-kurtosis matrix in vector form (collecting only the unique elements).}
#' \item{\code{Phi_shred}}{Partition on \code{Phi} (see reference).}
#' \item{\code{Psi_shred}}{Partition on \code{Psi} (see reference).}
#'
#'
#' @examples
#' 
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' X_moments <- estimate_sample_moments(X50[, 1:10])
#' 
#'
#' @importFrom PerformanceAnalytics M3.MM
#' @importFrom PerformanceAnalytics M4.MM
#' @importFrom stats cov
#' @export
estimate_sample_moments <- function(X, adjust_magnitude = FALSE) {
  M3.mat2vec <- get("M3.mat2vec", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)
  M3.vec2mat <- get("M3.vec2mat", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)
  M4.vec2mat <- get("M4.vec2mat", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)
  
  
  N <- ncol(X)
  mu  <- colMeans(X)
  Sgm <- cov(X)
  Phi <- M3.MM(X, as.mat = FALSE)
  Psi <- M4.MM(X, as.mat = FALSE)
  
  if (adjust_magnitude) {
    tmp_list <- list(mu = mu, Sgm = Sgm, Phi = Phi, Psi = Psi)
    attr(tmp_list, "type") <- "X_sample_moments"
    d <- abs(eval_portfolio_moments(w = rep(1/N, N), X_statistics = tmp_list))
    mu  <- mu  / d[1]
    Sgm <- Sgm / d[2]
    Phi <- Phi / d[3]
    Psi <- Psi / d[4]
  }
  
  # compute shred versions from mat versions
  Phi_mat <- M3.vec2mat(Phi, N)
  Phi_shred <- lapply(1:N, function(i) Phi_mat[, (1:N)+N*(i-1)])
  
  Psi_mat <- M4.vec2mat(Psi, N)
  Psi_shred <- list()
  for (i in 1:N) {
    tmp <- Psi_mat[, (1:N^2)+N^2*(i-1)]
    Psi_shred[[i]] <- M3.mat2vec(tmp)
    gc()
  }
  gc()
  
  
  list_to_return <- list(mu = mu, Sgm = Sgm, Phi_mat = Phi_mat, Psi_mat = Psi_mat, Phi = Phi, Psi = Psi, Phi_shred = Phi_shred, Psi_shred = Psi_shred)
  attr(list_to_return, "type") <- "X_sample_moments"
  return(list_to_return)
}




#' @title Evaluate first four moments of a given portfolio
#'
#' @description Evaluate first four moments of a given portfolio's return, namely,
#' mean, variance, skewness, and kurtosis.
#'
#' @author Rui Zhou, Xiwen Wang, and Daniel P. Palomar
#'
#' @references
#' R. Zhou and D. P. Palomar, "Solving High-Order Portfolios via Successive Convex Approximation Algorithms," 
#' in \emph{IEEE Transactions on Signal Processing}, vol. 69, pp. 892-904, 2021.
#' <doi:10.1109/TSP.2021.3051369>.
#' 
#' X. Wang, R. Zhou, J. Ying, and D. P. Palomar, "Efficient and Scalable High-Order Portfolios Design via Parametric Skew-t Distribution," 
#' Available in arXiv, 2022. <https://arxiv.org/pdf/2206.02412v1.pdf>.
#'
#' @param w Numerical vector with portfolio weights.
#' @param var Argument characterizing the constituents assets. 
#'            Either the sample parameters as obtained by function \code{\link{estimate_sample_moments}()} or
#'            the multivariate skew t parameters as obtained by function \code{\link{estimate_skew_t}()}.
#' 
#' @return Four moments of the given portfolio.
#' 
#' @examples
#' 
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' # nonparametric case
#' X_moments <- estimate_sample_moments(X50[, 1:10])
#' w_moments <- eval_portfolio_moments(w = rep(1/10, 10), X_statistics = X_moments)
#' 
#' # parametric case (based on the multivariate skew t distribution)
#' X_skew_t_params <- estimate_skew_t(X50[, 1:10])
#' w_moments <- eval_portfolio_moments(w = rep(1/10, 10), X_statistics = X_skew_t_params)
#' 
#'
#' @export
eval_portfolio_moments <- function(w, X_statistics) {
  
  if (attr(X_statistics, "type") == "X_sample_moments") {
    portm3 <- get("portm3", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
    portm4 <- get("portm4", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
    
    return(c(mean     = as.numeric(w %*% X_statistics$mu),
             variance = as.numeric(w %*% X_statistics$Sgm %*% w),
             skewness = as.numeric(portm3(w = w, M3 = X_statistics$Phi)),
             kurtosis = as.numeric(portm4(w = w, M4 = X_statistics$Psi))))
  } else if (attr(X_statistics, "type") == "X_skew_t_params") {
    # Xiwen
  } else
    stop("Unknown type of argument ", dQuote(X_statistics), " : it should be returned from either", dQuote("estimate_sample_moments()"), 
    " or ", dQuote("estimate_skew_t()"))
}


