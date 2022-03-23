# This is a library of auxiliary functions, which might be useful for other exported functions of this package.
# They should be invisible to package users

# looking for the closest positive semidefinite approximation matrix for a given hessain matrix --------------------------
# .apprxHessian <- function(hsn) {
#   eig_decomp <- eigen(hsn)
#   eig_decomp$vectors %*% diag(pmax(eig_decomp$values, 0)) %*% t(eig_decomp$vectors)
# }

.apprxHessian <- function(hsn, decomp = FALSE) {
  eig_decomp <- eigen(hsn)
  n <- max(1, sum(eig_decomp$values > 0))
  L <- eig_decomp$vectors[, 1:n, drop = FALSE] %*% diag(sqrt(pmax(eig_decomp$values[1:n], 0)), n)
  if (decomp) 
    return(list("hsn" = L%*%t(L), "L" = L))
  else
    return(L%*%t(L))
}

# upper bound for eigenvalue of Hessian of skewness (when leverage == 1) -------------------------------------------------
#' @importFrom magrittr %>% multiply_by
.maxEigHsnS <- function(S, N, func = "max") {
  M3.vec2mat <- get("M3.vec2mat", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
  
  if (is.vector(S)) S <- M3.vec2mat(S, N)
  S <- abs(S)
  if (func == "max")
    res <- do.call(pmax, lapply(1:N, function(i) S[, .idx_mask(i, N)])) %>% rowSums() %>% max() %>% multiply_by(6)
  if (func == "sum")
    res <- Reduce("+", lapply(1:N, function(i) S[, .idx_mask(i, N)])) %>% rowSums() %>% max() %>% multiply_by(6) 
  return(res)
}

#' @importFrom magrittr %>%
#  upper bound for eigenvalue of Hessian of kurtosis (when leverage == 1) ------------------------------------------------
.maxEigHsnK <- function(K, N, func = "max") {
  M4.vec2mat <- get("M4.vec2mat", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
  
  if (is.vector(K)) K <- M4.vec2mat(K, N)
  K <- abs(K)
  if (func == "max")
    res <- do.call(pmax, lapply(1:(N^2), function(i) K[, .idx_mask(i, N)])) %>% rowSums() %>% max() %>% multiply_by(12)
  if (func == "sum")
    res <- Reduce("+", lapply(1:(N^2), function(i) K[, .idx_mask(i, N)])) %>% rowSums() %>% max() %>% multiply_by(12)
  return(res)
}

# pruducing the mask -----------------------------------------------------------------------------------------------------
.idx_mask <- function(i, N) 1:N + (i - 1)*N


# formating a quadratic inequality to cope with SOCP solver interface ----------------------------------------------------
.QCQP2SOCP <- function(Q, q, l, L = NULL) {
  # xQx + qx + l <= 0
  if (is.null(L)) {
    tmp <- eigen(Q)
    L <- tmp$vectors %*% diag(sqrt(tmp$values))
  }
  N <- ncol(L)
  list(
    "G" = rbind(q/2, t(L), q/2),
    "h" = c((1-l)/2, rep(0, N), -(1+l)/2)
  )
}

# a benchmark for MVSK portfolio, implemented with package nloptr --------------------------------------------------------
#' @import nloptr
.MVSKnloptr <- function(lmd = rep(1, 4), mom_params, w0 = rep(1/length(mom_params$mu), length(mom_params$mu)),
                        stopval = -Inf, xtol_rel = 1e-6, ftol_rel = 1e-6) {
  portm3 <- get("portm3", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
  portm4 <- get("portm4", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
  derportm3 <- get("derportm3", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
  derportm4 <- get("derportm4", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)  
  
  # extract moment parameters
  mu  <- mom_params$mu
  Sgm <- mom_params$Sgm
  Phi <- mom_params$Phi
  Psi <- mom_params$Psi
  
  N <- length(w0)
  f <- function(w) {
    - lmd[1] * as.numeric(w %*% mu) +
      + lmd[2] * as.numeric(w %*% Sgm %*% w) +
      - lmd[3] * as.numeric(portm3(w = w, M3 = Phi)) +
      + lmd[4] * as.numeric(portm4(w = w, M4 = Psi))
  }
  
  grad_f <- function(w) {
    - lmd[1] * mu + 
      + lmd[2] * as.numeric(2 * Sgm %*% w) +
      - lmd[3] * derportm3(w, Phi) +
      + lmd[4] * derportm4(w, Psi) 
  }
  
  constraints <- function(w) {
    list("constraints" = sum(w) - 1,
         "jacobian" = rep(1, N))
  }
  local_opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                      "xtol_rel" = xtol_rel,
                      "ftol_rel" = ftol_rel)
  
  start_time <- proc.time()[3]
  res0 <- nloptr(x0 = rep(1/N, N),
                 eval_f = f,
                 eval_grad_f = grad_f,
                 lb = rep(0, N),
                 ub = rep(1, N),
                 eval_g_eq = constraints,
                 opts = list("algorithm" = "NLOPT_LD_AUGLAG",
                             "local_opts" = local_opts,
                             "xtol_rel" = xtol_rel,
                             "ftol_rel" = ftol_rel,
                             "stopval" = stopval,
                             "maxeval" = 1e4,
                             "print_level" = 0,
                             "check_derivatives" = FALSE,
                             "check_derivatives_print" = "all"))
  tm <- as.numeric(proc.time()[3] - start_time)
  # browser()
  return(list(
    "w" = res0$solution,
    "obj" = f(res0$solution),
    "time" = tm
  ))
  
}
