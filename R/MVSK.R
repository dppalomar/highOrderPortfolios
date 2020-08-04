#' @title Design high-order portfolio based on weighted linear combination of first four moments
#' 
#' @description Design high-order portfolio based on weighted linear combination of first four moments
#' (i.e., mean, variance, skewness, and kurtosis):
#' \preformatted{
#'   minimize     - lmd1*(w'*mu) + lmd2*(w'*Sigma*w) 
#'                - lmd3*(w'*Phi*w*w) + lmd4*(w'*Psi*w*w*w)
#'   subject to   ||w||_1 <= leverage, sum(w) == 1,
#' }
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' Rui Zhou and and Daniel P. Palomar, “Solving High-Order Portfolios via Successive Convex Approximation Algorithms,”
#' under review, 2020. <https://arxiv.org/abs/2008.00863>
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
#' \item{\code{cpu_time}}{Time usage with iteration.}
#' \item{\code{objs}}{Function objective with iteration.}
#' \item{\code{convergence}}{Boolean flag to indicate whether or not the optimization converged.}
#' \item{\code{moments}}{Moments of portfolio return at optimal portfolio weights.}
#'
#' @examples
#' \donttest{
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' # estimate moments
#' X_moments <- estimate_moments(X50)
#' 
#' # decide moment weights
#' xi <- 10
#' lmd <- c(1, xi/2, xi*(xi+1)/6, xi*(xi+1)*(xi+2)/24)
#' 
#' # portfolio optimization
#' sol <- design_MVSK_portfolio(lmd, X_moments)
#' }
#' 
#' @import PerformanceAnalytics
#' @importFrom utils tail
#' @export
design_MVSK_portfolio <- function(lmd = rep(1, 4), X_moments, 
                                  w_init = rep(1/length(X_moments$mu), length(X_moments$mu)), 
                                  leverage = 1, method = c("Q-MVSK", "MM", "DC"),
                                  tau_w = 0, gamma = 1, zeta = 1e-8, maxiter = 1e2, ftol = 1e-5, wtol = 1e-4, stopval = -Inf) {
  # error control
  if (leverage < 1) stop("leverage must be no less than 1.")
  if (leverage != 1) stop("Support for leverage > 1 is coming in next version.")
  
  # argument handling
  method <- match.arg(method)
  if (method != "Q-MVSK") {gamma = 1; zeta = 0}
  
  # extract moment parameters
  mu  <- X_moments$mu
  Sgm <- X_moments$Sgm
  Phi <- X_moments$Phi
  Psi <- X_moments$Psi
  Phi_shred <- X_moments$Phi_shred
  Psi_shred <- X_moments$Psi_shred
  
  # browser()
  N <- length(mu)
  Amat <- t(rbind(matrix(1, 1, N), diag(N)))
  bvec <- cbind(c(1, rep(0, N)))
  w <- w_init
  
  cpu_time <- c(0)
  objs  <- c() 
  getgrad <- function(w) rbind(mu, 2*w%*%Sgm, as.vector(PerformanceAnalytics:::derportm3(w, Phi)), as.vector(PerformanceAnalytics:::derportm4(w, Psi)))  # gradients computing function
  obj <- function() sum(lmd * as.vector(grads %*% w) / c(-1, 2, -3, 4))  # objective computing function, grads must be prepared before
  
  # when method is "MM" or "DC"
  if (method == "MM") rho <- leverage*lmd[3]*.maxEigHsnS(S = Phi, N = N, func = "max") + leverage^2*lmd[4]*.maxEigHsnK(K = Psi, N = N, func = "max")
  if (method == "DC") rho <- leverage*lmd[3]*.maxEigHsnS(S = Phi, N = N, func = "sum") + leverage^2*lmd[4]*.maxEigHsnK(K = Psi, N = N, func = "sum")
  
  # browser()
  start_time <- proc.time()[3]
  
  # compute current gradient and objective
  H3 <- 6 * sapply(Phi_shred, function(x) x%*%w)
  H4 <- 4 * sapply(Psi_shred, function(x) PerformanceAnalytics:::derportm3(w, x))
  H34 <- - lmd[3] * H3 + lmd[4] * H4
  grads <- rbind(mu, 2*w%*%Sgm, w%*%H3/2, w%*%H4/3)
  objs <- c(objs, obj())
  
  # SCA outer loop
  for (iter in 1:maxiter) {
    
    # record previous w
    w_old <- w
    
    ## construct QP approximation problem (the symbol and scale is adjusted to match the format of solver quadprog::solve.QP)
    if (method == "Q-MVSK") {  # method == "Q-MVSK
      H_ncvx <- .apprxHessian(H34)
      Qk <- 2*lmd[2]*Sgm + H_ncvx + diag(tau_w, N)
      qk <- lmd[1]*mu + lmd[3]*grads[3, ] - lmd[4]*grads[4, ] + H_ncvx%*%w + tau_w*w
    }
    if (method == "MM") {  # method == "MM"
      Qk <- 2*lmd[2]*Sgm + diag(rho, N)
      qk <- lmd[1]*mu + lmd[3]*grads[3, ] - lmd[4]*grads[4, ] + rho*w
    }
    if (method == "DC") {  # method == "DC"
      Qk <- diag(rho, N)
      qk <- rho*w + lmd[1]*grads[1, ] - lmd[2]*grads[2, ] + lmd[3]*grads[3, ] - lmd[4]*grads[4, ]
    }
    
    
    
    # browser()
    # solve the QP problem
    w_hat <- quadprog::solve.QP(Dmat = Qk, dvec = cbind(qk), Amat = Amat, bvec = bvec, meq = 1)$solution
    
    # update w
    w <- w + gamma * (w_hat - w)
    gamma <- gamma * (1 - zeta * gamma)
    
    # recording ...
    cpu_time <- c(cpu_time, proc.time()[3] - start_time) 
    
    # Hessian matrix
    H3 <- 6 * sapply(Phi_shred, function(x) x%*%w)
    H4 <- 4 * sapply(Psi_shred, function(x) PerformanceAnalytics:::derportm3(w, x))
    H34 <- - lmd[3] * H3 + lmd[4] * H4
    
    # recovery gradients from Hessian information
    grads <- rbind(mu, 2*w%*%Sgm, w%*%H3/2, w%*%H4/3)
    
    # record current objective
    objs <- c(objs, obj())
    
    # judge convergence
    # has_w_converged <- all(abs(w - w_old) <= .5 * wtol )
    # has_w_converged <- norm(w - w_old, "2") <= wtol * norm(w_old, "2")
    has_w_converged <- all(abs(w - w_old) <= .5 * wtol * (abs(w) + abs(w_old)))
    has_f_converged <- abs(diff(tail(objs, 2))) <= .5 * ftol * sum(abs(tail(objs, 2)))
    has_cross_stopval <- tail(objs, 1) <= stopval
    
    if (has_w_converged || has_f_converged || has_cross_stopval) break
  }
  
  return(list(
    "w"           = w,
    "cpu_time"    = cpu_time,
    "objs"        = objs,
    "convergence" = !(iter == maxiter),
    "moments"     = as.vector(grads %*% w) / c(1, 2, 3, 4)
  ))
  
}

