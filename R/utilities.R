#' @title Estimate the parameters of skew-t distribution from multivariate observations
#'
#' @description Using the package fitHeavyTail to estimate the parameters of ghMST distribution from multivariate observations, namely,
#' location vector (mu), skewness vector (gamma), scatter matrix (scatter), degree of freedom (nu), parameters a,
#' and the Cholesky decomposition of the scatter matrix (chol_Sigma).
#'
#' @author Xiwen Wang, Rui Zhou, and Daniel P. Palomar
#'
#' @references
#' Aas, Kjersti and Ingrid Hobæk Haff. "The generalized hyperbolic skew student’st-distribution,"
#' Journal of financial econometrics, pp. 275-309, 2006.
#' 
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param initial List of initial values of the parameters for the iterative estimation method.
#'                Possible elements include:
#'                \itemize{\item{\code{nu}: default is \code{4},}
#'                         \item{\code{mu}: default is the data sample mean,}
#'                         \item{\code{gamma}: default is the sample skewness vector,}
#'                         \item{\code{scatter}: default follows from the scaled sample covariance matrix,}
#'                         }
#' @param nu_lb Minimum value for the degree of freedom to maintain the existence of high-order moments (default is \code{9}).
#' @param max_iter Integer indicating the maximum number of iterations for the iterative estimation
#'                 method (default is \code{100}).
#' @param ptol Positive number indicating the relative tolerance for the change of the variables
#'             to determine convergence of the iterative method (default is \code{1e-3}).
#' @param ftol Positive number indicating the relative tolerance for the change of the log-likelihood
#'             value to determine convergence of the iterative method (default is \code{Inf}, so it is
#'             not active). Note that using this argument might have a computational cost as a convergence
#'             criterion due to the computation of the log-likelihood (especially when \code{X} is high-dimensional).
#' @param PXEM Logical value indicating whether to use the parameter expansion (PX) EM method to accelerating the convergence.
#' @param return_iterates Logical value indicating whether to record the values of the parameters (and possibly the
#'                        log-likelihood if \code{ftol < Inf}) at each iteration (default is \code{FALSE}).
#' @param verbose Logical value indicating whether to allow the function to print messages (default is \code{FALSE}).
#'
#' @return A list containing the following elements:
#'         \item{\code{mu}}{Location vector estimate (not the mean).}
#'         \item{\code{gamma}}{Skewness vector estimate.}
#'         \item{\code{scatter}}{Scatter matrix estimate.}
#'         \item{\code{nu}}{Degrees of freedom estimate.}
#'         \item{\code{chol_Sigma}}{Choleski decomposition of the Scatter matrix estimate.}
#'         \item{\code{a}}{A list of coefficients useful for later computation}
#'
#' @examples
#' library(highOrderPortfolios)
#' data("X50")
#' X_skew_t_params <- estimate_skew_t(X50)
#'
#' @import fitHeavyTail
#' @export
estimate_skew_t <- function(X, initial = NULL, nu_lb = 9, max_iter = 100, ptol = 1e-3, ftol = Inf,
                            PXEM = TRUE, return_iterates = FALSE, verbose = FALSE) { 
  options(nu_min = nu_lb)
  X_skew_t_params <- fitHeavyTail::fit_mvst(X, max_iter = max_iter, ptol = ptol, ftol = ftol,
                                            PXEM = PXEM, return_iterates = return_iterates, verbose = verbose)
  return_paramters <- list()
  return_paramters$mu    <- X_skew_t_params$mu
  return_paramters$nu    <- X_skew_t_params$nu
  return_paramters$gamma <- X_skew_t_params$gamma
  return_paramters$scatter <- X_skew_t_params$scatter
  return_paramters$chol_Sigma <- chol(return_paramters$scatter)
  
  ## Compute a given paramters of skew-t model
  compute_a <- function(X_skew_t_params) {
    a11 <- (X_skew_t_params$nu)/(X_skew_t_params$nu-2)
    a21 <- (X_skew_t_params$nu)/(X_skew_t_params$nu-2)
    a22 <- (2 * X_skew_t_params$nu ** 2)/((X_skew_t_params$nu - 2) ** 2 * (X_skew_t_params$nu - 4))
    a31 <- 16 * (X_skew_t_params$nu ** 3)/(((X_skew_t_params$nu - 2) ** 3) * (X_skew_t_params$nu - 4) * (X_skew_t_params$nu - 6))
    a32 <- 6 * (X_skew_t_params$nu ** 2)/(((X_skew_t_params$nu -2) ** 2) * (X_skew_t_params$nu - 4))
    a41 <- (12 * X_skew_t_params$nu + 120) * (X_skew_t_params$nu ** 4)/ (((X_skew_t_params$nu - 2) ** 4) * (X_skew_t_params$nu - 4) * (X_skew_t_params$nu - 6) * (X_skew_t_params$nu - 8))
    a42 <- (2 * X_skew_t_params$nu + 4) * (X_skew_t_params$nu ** 3)/(((X_skew_t_params$nu - 2) ** 3) * (X_skew_t_params$nu - 4) * (X_skew_t_params$nu - 6)) * 6
    a43 <- (X_skew_t_params$nu ** 2)/((X_skew_t_params$nu - 2) * (X_skew_t_params$nu - 4)) * 3
    return(list("a11" = a11, "a21" = a21, "a22" = a22,
                "a31" = a31, "a32" = a32,
                "a41" = a41, "a42" = a42, "a43" = a43))
  }
  
  return_paramters$a <- compute_a(return_paramters)
  attr(return_paramters, "type") <- "X_skew_t_params"
  return(return_paramters)
}




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
#'                         Note: this is specially designed for the function \code{\link{design_MVSKtilting_portfolio_via_sample_moments}()}.
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
#' @param X_statistics Argument characterizing the constituents assets. 
#'                     Either the sample parameters as obtained by function \code{\link{estimate_sample_moments}()} or
#'                     the multivariate skew t parameters as obtained by function \code{\link{estimate_skew_t}()}.
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
    wT_mu <- t(w) %*% X_statistics$mu
    wT_gamma <- t(w) %*% X_statistics$gamma
    wT_Sigma_w <- sum((X_statistics$chol_Sigma %*% w) ** 2)
    return(c(mean     = as.numeric(wT_mu + X_statistics$a$a11 * wT_gamma),
             variance = as.numeric(X_statistics$a$a21 * wT_Sigma_w  + X_statistics$a$a22 * ((wT_gamma)**2) ),
             skewness = as.numeric(X_statistics$a$a31 * ((wT_gamma) ** 3) + X_statistics$a$a32 * (wT_gamma) * ( wT_Sigma_w) ),
             kurtosis = as.numeric(X_statistics$a$a41 * ((wT_gamma) ** 4) + X_statistics$a$a42 * (wT_Sigma_w) * ((wT_gamma) ** 2) + X_statistics$a$a43 * ((wT_Sigma_w) ** 2))))
  } else
    stop("Unknown type of argument ", dQuote(X_statistics), " : it should be returned from either", dQuote("estimate_sample_moments()"), 
    " or ", dQuote("estimate_skew_t()"))
}



