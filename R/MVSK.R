#' @title Design high-order portfolio based on weighted linear combination of first four moments
#' 
#' @description Design high-order portfolio based on weighted linear combination of first four moments
#' (i.e., mean, variance, skewness, and kurtosis):
#' \preformatted{
#'   minimize     - lmd1*(w'*mu) + lmd2*(w'*Sigma*w) 
#'                - lmd3*(w'*Phi*w*w) + lmd4*(w'*Psi*w*w*w)
#'   subject to   ||w||_1 <= leverage, sum(w) == 1.
#' }
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' R. Zhou and D. P. Palomar, "Solving High-Order Portfolios via Successive Convex Approximation Algorithms," 
#' in \emph{IEEE Transactions on Signal Processing}, vol. 69, pp. 892-904, 2021.
#' <doi:10.1109/TSP.2021.3051369>.
#' 
#' X. Wang, R. Zhou, J. Ying, and D. P. Palomar, "Efficient and Scalable High-Order Portfolios Design via Parametric Skew-t Distribution," 
#' Available in arXiv, 2022. <https://arxiv.org/pdf/2206.02412v1.pdf>.
#' 
#' @param lmd Numerical vector of length 4 indicating the weights of first four moments.
#' @param X_moments List of moment parameters, see \code{\link{estimate_moments}}.
#' @param w_init Numerical vector indicating the initial value of portfolio weights.
#' @param leverage Number (>= 1) indicating the leverage of portfolio.
#' @param method String indicating the algorithm method, must be one of: "Q-MVSK", "MM", "DC".
#' @param tau_w Number (>= 0) guaranteeing the strong convexity of approximating function.
#' @param gamma Number (0 < gamma <= 1) indicating the initial value of gamma.
#' @param zeta Number (0 < zeta < 1) indicating the diminishing paramater of gamma.
#' @param maxiter Positive integer setting the maximum iteration.
#' @param ftol Positive number setting the convergence criterion of function objective.
#' @param wtol Positive number setting the convergence criterion of portfolio weights.
#' @param stopval Number setting the stop value of objective.
#' 
#' @return A list containing the following elements:
#' \item{\code{w}}{Optimal portfolio vector.}
#' \item{\code{cpu_time_vs_iterations}}{Time usage over iterations.}
#' \item{\code{objfun_vs_iterations}}{Objective function over iterations.}
#' \item{\code{iterations}}{Iterations index.}
#' \item{\code{convergence}}{Boolean flag to indicate whether or not the optimization converged.}
#' \item{\code{moments}}{Moments of portfolio return at optimal portfolio weights.}
#'
#' @examples
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' # estimate moments
#' X_moments <- estimate_sample_moments(X50[, 1:10])
#' 
#' # decide moment weights
#' xi <- 10
#' lmd <- c(1, xi/2, xi*(xi+1)/6, xi*(xi+1)*(xi+2)/24)
#' 
#' # portfolio optimization
#' sol <- design_MVSK_portfolio_via_sample_moments(lmd, X_moments)
#' 
#' @importFrom utils tail
#' @import PerformanceAnalytics
#' @export
design_MVSK_portfolio_via_sample_moments <- function(lmd = rep(1, 4), X_moments, 
                                                     w_init = rep(1/length(X_moments$mu), length(X_moments$mu)), 
                                                     leverage = 1, method = c("Q-MVSK", "MM", "DC"),
                                                     tau_w = 0, gamma = 1, zeta = 1e-8, maxiter = 1e2, ftol = 1e-5, wtol = 1e-4, stopval = -Inf) {
  method <- match.arg(method)
  derportm3 <- get("derportm3", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
  
  # error control
  if (attr(X_moments, "type") != "X_sample_moments")
    stop("Argument X_moments is not of type ", dQuote("X_sample_moments"), ". It should be returned from function ", dQuote("estimate_sample_moments()"), ".")
  if (leverage < 1) stop("leverage must be no less than 1.")
  if (leverage != 1) stop("Support for leverage > 1 is coming in next version.")
  
  # prep
  N <- length(X_moments$mu)
  Amat <- t(rbind(matrix(1, 1, N), diag(N)))
  bvec <- cbind(c(1, rep(0, N)))
  
  fun_eval <- function() {
    return_list <- list()

    if (method == "Q-MVSK") {
      return_list$H3 <- 6 * sapply(X_moments$Phi_shred, function(x) x%*%w)
      return_list$H4 <- 4 * sapply(X_moments$Psi_shred, function(x) derportm3(w, x))
      return_list$H34 <- - lmd[3] * return_list$H3 + lmd[4] * return_list$H4
      return_list$jac <- rbind("grad1" = X_moments$mu, "grad2" = 2 * c(X_moments$Sgm %*% w), "grad3" = (1/2) * c(return_list$H3 %*% w), "grad4" = (1/3) * c(return_list$H4 %*% w))
      
    } else {
      w_kron_w <- kronecker(w, w)
      return_list$jac <- rbind("grad1" = X_moments$mu, "grad2" = 2 * c(X_moments$Sgm %*% w), "grad3" = 3 * c(X_moments$Phi_mat %*% w_kron_w), "grad4" = 4 * c(X_moments$Psi_mat %*% kronecker(w_kron_w, w)))
    }
    return_list$obj <- sum(lmd * as.vector(return_list$jac %*% w) / c(-1, 2, -3, 4))
    
    return(return_list)
  }
  
  # initialization
  start_time <- proc.time()[3]
  if (method == "MM") {
    gamma = 1; zeta = 0
    rho <- leverage*lmd[3]*.maxEigHsnS(S = X_moments$Phi, N = N, func = "max") + leverage^2*lmd[4]*.maxEigHsnK(K = X_moments$Psi, N = N, func = "max")
  }
  if (method == "DC") {
    gamma = 1; zeta = 0
    rho <- leverage*lmd[3]*.maxEigHsnS(S = X_moments$Phi, N = N, func = "sum") + leverage^2*lmd[4]*.maxEigHsnK(K = X_moments$Psi, N = N, func = "sum")
  }
  w <- w_init
  cpu_time <- c(0)
  objs  <- c()  
  fun_k <- fun_eval()
  objs <- c(objs, fun_k$obj)

  #
  # SCA outer loop
  #
  for (iter in 1:maxiter) {
    # record previous w
    w_old <- w
    
    ## construct QP approximation problem (the symbol and scale is adjusted to match the format of solver quadprog::solve.QP)
    switch(method, 
           "Q-MVSK" = {
             H_ncvx <- .apprxHessian(fun_k$H34)
             Qk <- 2*lmd[2]*X_moments$Sgm + H_ncvx + diag(tau_w, N)
             qk <- lmd[1]*X_moments$mu + lmd[3]*fun_k$jac[3, ] - lmd[4]*fun_k$jac[4, ] + H_ncvx%*%w + tau_w*w
           },
           "MM" = {
             Qk <- 2*lmd[2]*X_moments$Sgm + diag(rho, N)
             qk <- lmd[1]*X_moments$mu + lmd[3]*fun_k$jac[3, ] - lmd[4]*fun_k$jac[4, ] + rho*w
           },
           "DC" = {
             Qk <- diag(rho, N)
             qk <- rho*w + lmd[1]*fun_k$jac[1, ] - lmd[2]*fun_k$jac[2, ] + lmd[3]*fun_k$jac[3, ] - lmd[4]*fun_k$jac[4, ]
           },
           stop("Method unknown")
          )

    # solve the QP problem
    w_hat <- quadprog::solve.QP(Dmat = Qk, dvec = cbind(qk), Amat = Amat, bvec = bvec, meq = 1)$solution
    
    # update w
    w <- w + gamma * (w_hat - w)
    gamma <- gamma * (1 - zeta * gamma)
    
    # recording...
    cpu_time <- c(cpu_time, proc.time()[3] - start_time) 
    fun_k <- fun_eval()
    objs <- c(objs, fun_k$obj)
    
    # termination criterion
    # has_w_converged <- all(abs(w - w_old) <= .5 * wtol )
    # has_w_converged <- norm(w - w_old, "2") <= wtol * norm(w_old, "2")
    has_w_converged <- all(abs(w - w_old) <= .5 * wtol * (abs(w) + abs(w_old)))
    has_f_converged <- abs(diff(tail(objs, 2))) <= .5 * ftol * sum(abs(tail(objs, 2)))
    has_cross_stopval <- tail(objs, 1) <= stopval
    
    if (has_w_converged || has_f_converged || has_cross_stopval) break
  }
  
  return(list(
    "w"                      = w,
    "cpu_time_vs_iterations" = cpu_time,
    "objfun_vs_iterations"   = objs,
    "iterations"             = 0:iter,
    "convergence"            = !(iter == maxiter),
    "moments"                = as.vector(fun_k$jac %*% w) / c(1, 2, 3, 4)
  ))
}

