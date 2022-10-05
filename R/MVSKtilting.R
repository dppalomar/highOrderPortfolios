#' @title Design high-order portfolio by tilting a given portfolio to the MVSK efficient frontier
#' 
#' @description Design high-order portfolio by tilting a given portfolio to the MVSK efficient frontier
#' (i.e., mean, variance, skewness, and kurtosis):
#' \preformatted{
#'   minimize     - delta
#'                m1(w) >= m1(w0) + delta*d1
#'                m2(w) <= m2(w0) - delta*d2
#'                m3(w) >= m3(w0) + delta*d3
#'                m4(w) <= m4(w0) - delta*d4
#'                (w-w0)'Sigma(w-w0) <= kappa^2
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
#' @inheritParams design_MVSK_portfolio_via_sample_moments
#' @param d Numerical vector of length 4 indicating the weights of first four moments.
#' @param w0 Numerical vector indicating the reference portfolio vector.
#' @param w0_moments Numerical vector indicating the reference moments. 
#' @param kappa Number indicating the maximum tracking error volatility. 
#' @param tau_delta Number (>= 0) guaranteeing the strong convexity of approximating function.
#' @param theta Number (0 < theta < 1) setting the combination coefficient when enlarge feasible set.
#' 
#' @return A list containing the following elements:
#' \item{\code{w}}{Optimal portfolio vector.}
#' \item{\code{delta}}{Maximum tilting distance of the optimal portfolio.}
#' \item{\code{cpu_time_vs_iterations}}{Time usage over iterations.}
#' \item{\code{objfun_vs_iterations}}{Objective function over iterations.}
#' \item{\code{iterations}}{Iterations index.}
#' \item{\code{moments}}{Moments of portfolio return at optimal portfolio weights.}
#' \item{\code{improvement}}{The relative improvement of moments of designed portfolio w.r.t. the reference portfolio.}
#'
#' @examples
#' 
#' library(highOrderPortfolios)
#' data(X50)
#' 
#' # estimate moments
#' X_moments <- estimate_sample_moments(X50[, 1:10])
#' 
#' # decide problem setting
#' w0 <- rep(1/10, 10)
#' w0_moments <- eval_portfolio_moments(w0, X_moments)
#' d <- abs(w0_moments) 
#' kappa <- 0.3 * sqrt(w0 %*% X_moments$Sgm %*% w0)
#' 
#' # portfolio optimization
#' sol <- design_MVSKtilting_portfolio_via_sample_moments(d, X_moments, w_init = w0, w0 = w0, 
#'                                                        w0_moments = w0_moments, kappa = kappa)
#' 
#' 
#' @importFrom utils tail
#' @import quadprog
#' @import lpSolveAPI
#' @import ECOSolveR
#' @import PerformanceAnalytics
#' @export
design_MVSKtilting_portfolio_via_sample_moments <- function(d = rep(1, 4), X_moments, 
                                         w_init = rep(1/length(X_moments$mu), length(X_moments$mu)), 
                                         w0 = w_init, w0_moments = NULL, 
                                         leverage = 1, kappa = 0, method = c("Q-MVSKT", "L-MVSKT"),
                                         tau_w = 1e-5, tau_delta = 1e-5, gamma = 1, zeta = 1e-8, maxiter = 1e2, ftol = 1e-5, wtol = 1e-5, 
                                         theta = 0.5, stopval = -Inf) {
  method <- match.arg(method)
  derportm3 <- get("derportm3", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)
  derportm4 <- get("derportm4", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)
  
  # error control
  # error control
  if (attr(X_moments, "type") != "X_sample_moments")
    stop("Argument X_moments is not of type ", dQuote("X_sample_moments"), ". It should be returned from function ", dQuote("estimate_sample_moments()"), ".")
  if (leverage < 1) stop("leverage must be no less than 1.")
  if (leverage != 1) stop("Support for leverage > 1 is coming in next version.")
  
  # prep
  N <- length(X_moments$mu)
  
  fun_eval <- function() {
    return_list <- list()
    
    return_list$H3 <- 6 * sapply(X_moments$Phi_shred, function(x) x%*%w)
    return_list$H4 <- 4 * sapply(X_moments$Psi_shred, function(x) derportm3(w, x))
    return_list$jac <- rbind("grad1" = X_moments$mu, "grad2" = 2 * c(X_moments$Sgm %*% w), "grad3" = (1/2) * c(return_list$H3 %*% w), "grad4" = (1/3) * c(return_list$H4 %*% w))
    return_list$w_moments <- as.vector(return_list$jac %*% w) / c(1, 2, 3, 4)
    return_list$obj <- -min((return_list$w_moments - w0_moments) / d * c(1, -1, 1, -1))
    
    return(return_list)
  }
  
  # initialization
  start_time <- proc.time()[3]
  w <- w_init
  delta <- 0
  if (is.null(w0_moments)) 
    w0_moments <- eval_portfolio_moments(w = w, X_statistics = X_moments)
  cpu_time <- c(0)
  objs  <- c()
  fun_k <- fun_eval()
  objs <- c(objs, fun_k$obj)

  if (method == "Q-MVSKT") {
    # options_ecos <- ECOSolveR::ecos.control(feastol = ftol, reltol = ftol, abstol = ftol)
    options_ecos <- ECOSolveR::ecos.control()
    A_basic = rbind(c(rep(1, N), 0, 0))  # equality constraint: sum(w) == 1
    b_basic = 1
    G_basic <- cbind(-diag(N+1), 0)      # inequality constraint: w >= 0, delta >= 0, (and t >= 0)
    h_basic <- rep(0, N+1)
    L2 <- .apprxHessian(X_moments$Sgm, TRUE)$L
  }
  if (method == "L-MVSKT") {  # initialize QP solver and LP solver, setting some constant parameters
    # QP solver: pre-setted parameter
    A_mat <- cbind(c(rep(1, N), 0), diag(N+1))
    Dmat <- diag(c(rep(tau_w, N), tau_delta))
    
    # LP solver:
    lp <- make.lp(6, N+2) # 1 euqality constraint, 5 inequality constraint
    set.objfn(lprec = lp, obj = c(rep(0, N), 0, 1))
    set.row(lprec = lp, row = 1, xt = c(rep(1, N), 0, 0))
    set.constr.type(lprec = lp, types = c("=", rep(">=", 5)))
    # set.bounds(lprec = lp, upper = c(rep(beta_w, N), beta_delta, Inf))
  }
  
  #
  # SCA outer loop
  #
  for (iter in 1:maxiter) {
    # record previous w and delta
    w_old <- w; delta_old <- delta
    
    # compute eta for enlarging the feasible set of approximating problem
    if (method == "Q-MVSKT") {  
      # approximate and decompose Hessian matrix
      tmp3 <- .apprxHessian(-fun_k$H3, TRUE); tmp4 <- .apprxHessian(fun_k$H4, TRUE)
      L3 <- tmp3$L; L4 <- tmp4$L; H3_app <- tmp3$hsn; H4_app <- tmp4$hsn
      gk <- (fun_k$w_moments - w0_moments) * c(-1, 1, -1, 1) + delta * d  
      if (all(gk[3:4] <= 0)) {  # only g_3 and g_4 need approximation
        eta <- 0
      } else {
        # tracking error (g_5) constraint
        tmpRef <- .QCQP2SOCP(q = c(-2*as.vector(w0%*%X_moments$Sgm), 0, 0),
                             l = as.numeric(w0%*%X_moments$Sgm%*%w0) - kappa^2,
                             L = rbind(rbind(L2, 0), 0))
        GRef <- tmpRef$G
        hRef <- tmpRef$h
        
        # first moment (g_1) constraint
        h1 = - w0_moments[1]  
        G1 <- rbind(c(-X_moments$mu, d[1], 0))
        
        # second moment (g_2) constraint
        tmp2 <- .QCQP2SOCP(q = c(rep(0, N), d[2], 0),
                           l = - w0_moments[2],
                           L = rbind(rbind(L2, 0), 0))
        G2 <- tmp2$G
        h2 <- tmp2$h
        
        # third moment (g_3) constraint
        tmp3 <- .QCQP2SOCP(q = c(-fun_k$jac[3, ]-as.vector(w%*%H3_app), d[3], -1), 
                           l = gk[3] + sum(fun_k$jac[3, ]*w) - d[3]*delta +as.numeric( w%*%H3_app%*%w/2),
                           L = rbind(rbind(L3/sqrt(2), 0), 0))
        G3 <- tmp3$G
        h3 <- tmp3$h
        
        # fourth moment (g_4) constraint
        tmp4 <- .QCQP2SOCP(q = c(fun_k$jac[4, ]-as.vector(w%*%H4_app), d[4], -1), 
                           l = gk[4] - sum(fun_k$jac[4, ]*w) - d[4]*delta + as.numeric(w%*%H4_app%*%w/2),
                           L = rbind(rbind(L4/sqrt(2), 0), 0))
        G4 <- tmp4$G
        h4 <- tmp4$h
        
        # solve problem
        sol_socp <- ECOSolveR::ECOS_csolve(c = c(rep(0, N+1), 1), 
                                           G = rbind(G_basic, G1, G2, G3, G4, GRef), 
                                           h = c(h_basic, h1, h2, h3, h4, hRef), 
                                           dims = list(l = length(h_basic) + length(h1), q = c(length(h2), length(h3), length(h4), length(hRef)), e = 0L), 
                                           A = A_basic, b = b_basic, control = options_ecos)
        eta <- theta*max(sol_socp$x[N+2], 0) + (1 - theta) * max(gk[3:4])
      }
    }
    if (method == "L-MVSKT") {
      f.con <- cbind(rbind(-2 * as.numeric((w-w0)%*%X_moments$Sgm), fun_k$jac * c(1, -1, 1, -1)), c(0, -d))
      gk <- c((w-w0)%*%X_moments$Sgm%*%(w-w0) - kappa^2, (fun_k$w_moments - w0_moments) * c(-1, 1, -1, 1) + delta * d) 
      f.rhs <- gk + f.con%*%c(w, delta)
      if ( all(gk <= 0) ) {
        eta <- 0
      } else {
        for (i in 1:5) set.row(lprec = lp, row = i+1, xt = c(f.con[i, ], 1))
        set.rhs(lprec = lp, b = c(1, f.rhs))  # TODO: check this 
        # browser()
        solve(lp)
        eta <- theta*get.objective(lp) + (1-theta)*max(gk)
      }
    }
    
    

    # solve the approximating problem
    if (method == "Q-MVSKT") {
      # the objective
      tmpObj <- .QCQP2SOCP(q = c(-tau_w*w, -tau_delta*delta - 1, -1),
                           l = tau_delta*delta^2/2 + tau_w*sum(w^2)/2,
                           L = diag(c(rep(sqrt(tau_w/2), N), sqrt(tau_delta/2), 0)))
      GObj <- tmpObj$G
      hObj <- tmpObj$h
      
      # tracking error (g_5) constraint
      tmpRef <- .QCQP2SOCP(q = c(-2*as.vector(w0%*%X_moments$Sgm), 0, 0),
                           l = as.numeric(w0%*%X_moments$Sgm%*%w0) - kappa^2,
                           L = rbind(rbind(L2, 0), 0))
      GRef <- tmpRef$G
      hRef <- tmpRef$h
      
      # first moment (g_1) constraint
      h1 = - w0_moments[1]  
      G1 <- rbind(c(-X_moments$mu, d[1], 0))
      
      # second moment (g_2) constraint
      tmp2 <- .QCQP2SOCP(q = c(rep(0, N), d[2], 0),
                         l = - w0_moments[2],
                         L = rbind(rbind(L2, 0), 0))
      G2 <- tmp2$G
      h2 <- tmp2$h
      
      # third moment (g_3) constraint
      tmp3 <- .QCQP2SOCP(q = c(-fun_k$jac[3, ]-as.vector(w%*%H3_app), d[3], 0),
                         l = gk[3] + sum(fun_k$jac[3, ]*w) - d[3]*delta +as.numeric( w%*%H3_app%*%w/2) - eta,
                         L = rbind(rbind(L3/sqrt(2), 0), 0))
      G3 <- tmp3$G
      h3 <- tmp3$h
      
      # fourth moment (g_4) constraint
      tmp4 <- .QCQP2SOCP(q = c(fun_k$jac[4, ]-as.vector(w%*%H4_app), d[4], 0),
                         l = gk[4] - sum(fun_k$jac[4, ]*w) - d[4]*delta + as.numeric(w%*%H4_app%*%w/2) - eta,
                         L = rbind(rbind(L4/sqrt(2), 0), 0))
      G4 <- tmp4$G
      h4 <- tmp4$h
      
      # solve problem
      sol_socp <- ECOSolveR::ECOS_csolve(c = c(rep(0, N+1), 1),
                                         G = rbind(G_basic, G1, GObj, GRef, G2, G3, G4),
                                         h = c(h_basic, h1, hObj, hRef, h2, h3, h4),
                                         dims = list(l = length(h_basic) + length(h1), q = c(length(hObj), length(hRef), length(h2), length(h3), length(h4)), e = 0L),
                                         A = A_basic, b = b_basic, control = options_ecos)
      
      w_hat <- sol_socp$x[1:N]
      delta_hat <- sol_socp$x[N+1]
    }
    
    if (method == "L-MVSKT") {
      bvec <- c(1, rep(0, N+1), f.rhs - eta)
      dvec <- c(tau_w*w, tau_delta*delta + 1)
      tmp <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = cbind(A_mat, t(f.con)), bvec = bvec, meq = 1)$solution
      w_hat <- tmp[1:N]
      delta_hat <- tmp[N+1]
    }
    
    
    # update w
    w <- w + gamma * (w_hat - w)
    delta <- delta + gamma * (delta_hat - delta)
    gamma <- gamma * (1 - zeta * gamma)
    
    # recording...
    cpu_time <- c(cpu_time, proc.time()[3] - start_time) 
    fun_k <- fun_eval()
    objs <- c(objs, fun_k$obj)

    # termination criterion
    # has_w_converged <- all(abs(w - w_old) <= .5 * wtol )
    has_w_converged <- norm(w - w_old, "2") <= wtol * norm(w_old, "2")
    has_f_converged <- abs(diff(tail(objs, 2))) <= ftol * abs(tail(objs, 1))
    has_cross_stopval <- tail(objs, 1) <= stopval
    
    if (has_w_converged || has_f_converged || has_cross_stopval) break
  }
  
  return(list(
    "w"                      = w,
    "delta"                  = delta,
    "cpu_time_vs_iterations" = cpu_time,
    "objfun_vs_iterations"   = objs,
    "iterations"             = 0:iter,
    "convergence"            = !(iter == maxiter),
    "moments"                = fun_k$w_moments,
    "improvement"            = (fun_k$w_moments - w0_moments) / d * c(1, -1, 1, -1)
  ))
  
  browser()  # this is necessary to avoid errors with ECOSOlveR package...
}

