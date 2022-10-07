#' @title Design MVSK portfolio without shorting based on the parameters of generalized hyperbolic skew-t distribution
#'
#' @description Design MVSK portfolio without shorting based on the parameters of generalized hyperbolic skew-t distribution:
#' \preformatted{
#'   minimize     - lambda1*phi1(w) + lambda2*phi2(w)
#'                - lambda3*phi3(w) + lambda4*phi4(w)
#'   subject to   w>=0, sum(w) == 1.
#' }
#'
#' @author Xiwen Wang, Rui Zhou and Daniel P. Palomar
#'
#' @references
#' X. Wang, R. Zhou, J. Ying, and D. P. Palomar, "Efficient and Scalable High-Order Portfolios Design via Parametric Skew-t Distribution,"
#' Available in arXiv, 2022. <https://arxiv.org/pdf/2206.02412.pdf>.
#'
#' @param lambda Numerical vector of length 4 indicating the weights of first four moments.
#' @param X_skew_t_params List of fitted parameters, including location vector, skewness vector, scatter matrix, and the degree of freedom,
#'                        see \code{\link{estimate_skew_t}()}.
#' @param w_init Numerical vector indicating the initial value of portfolio weights.
#' @param method String indicating the algorithm method, must be one of: "L-MVSK", "DC", "Q-MVSK", "SQUAREM", "RFPA", "PGD".
#' @param tau_w Number (>= 0) guaranteeing the strong convexity of approximating function.
#' @param beta Number (0 < beta < 1) decreasing the step size of the projected gradient methods.
#' @param tau Number (tau > 0) hyper-parameters for the fixed-point acceleration.
#' @param gamma Number (0 < gamma <= 1) indicating the initial value of gamma for the Q-MVSK method.
#' @param zeta Number (0 < zeta < 1) indicating the diminishing parameter of gamma for the Q-MVSK method.
#' @param maxiter Positive integer setting the maximum iteration.
#' @param initial_eta Initial eta for projected gradient methods
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
#' # estimate skew t distribution
#' X_skew_t_params <- estimate_skew_t(X50)
#' 
#' # decide moment weights
#' xi <- 10
#' lambda <- c(1, 4, 10, 20)
#' 
#' # portfolio optimization
#' sol <- design_MVSK_portfolio_via_skew_t(lambda, X_skew_t_params, method = "RFPA", tau = 10)
#'
#' @importFrom utils tail
#' @import quadprog
#' @export
design_MVSK_portfolio_via_skew_t <- function(lambda, X_skew_t_params,
                                             w_init = rep(1/length(X_skew_t_params$mu), length(X_skew_t_params$mu)),
                                             method = c("L-MVSK", "DC", "Q-MVSK", "SQUAREM", "RFPA", "PGD"), gamma = 1, zeta = 1e-8,
                                             tau_w = 0, beta = 0.5, tau = 1e5, initial_eta = 5, maxiter = 1e3, ftol = 1e-6, wtol = 1e-6, stopval = -Inf)  {
  # error control
  if (attr(X_skew_t_params, "type") != "X_skew_t_params")
    stop("Unknown type of argument ", dQuote("X_skew_t_params"), " : it should be returned from ", dQuote("estimate_skew_t()"))

  # re-format the lambda
  lambda_max <- max(lambda)
  lambda <- lambda/lambda_max
  # argument handling
  method <- match.arg(method)

  # prep
  N <- length(X_skew_t_params$mu)
  Amat <- t(rbind(matrix(1, 1, N), diag(N)))
  bvec <- cbind(c(1, rep(0, N)))

  # get the obj values from jac
  getobjs_skewt <- function(w, jac) {
    sum(lambda * as.vector(jac %*% w) / c(-1, 2, -3, 4))
  }

  fun_eval <- function(w) {
    return_list <- list()

    if (method == "Q-MVSK") {
      return_list$H2 <- 2 * (X_skew_t_params$a$a21 * X_skew_t_params$scatter + X_skew_t_params$a$a22 * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma))
      return_list$H3 <- (X_skew_t_params$a$a31 * 6 * as.numeric(t(w) %*% X_skew_t_params$gamma) * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma) + X_skew_t_params$a$a32 * (2 * X_skew_t_params$gamma %*% t(w) %*% X_skew_t_params$scatter + 2 * X_skew_t_params$scatter %*% w %*% t(X_skew_t_params$gamma) + 2 * as.numeric(t(w) %*% X_skew_t_params$gamma) * X_skew_t_params$scatter ))
      return_list$H4 <- (12 * X_skew_t_params$a$a41 * (as.numeric(t(w) %*% X_skew_t_params$gamma)**2) * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma)  + 2 * X_skew_t_params$a$a42 * (2 * as.numeric(t(w) %*% X_skew_t_params$gamma) * X_skew_t_params$scatter %*% w %*% t(X_skew_t_params$gamma) + ((as.numeric(t(w) %*% X_skew_t_params$gamma))**2) * X_skew_t_params$scatter + 2 * as.numeric(t(w) %*% X_skew_t_params$gamma) * X_skew_t_params$gamma %*% t(w) %*% X_skew_t_params$scatter + as.numeric(t(w) %*% X_skew_t_params$scatter %*% w) * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma)    ) + 4 * X_skew_t_params$a$a43 * (as.numeric(t(w) %*% X_skew_t_params$scatter %*% w) * X_skew_t_params$scatter + 2 * X_skew_t_params$scatter %*% w %*% t(w) %*% X_skew_t_params$scatter))
      return_list$jac <- rbind(as.vector(X_skew_t_params$mu + X_skew_t_params$a$a11 * X_skew_t_params$gamma),
                        w%*%return_list$H2, w%*%return_list$H3/2, w%*%return_list$H4/3)
      return_list$obj <- sum(lambda * as.vector(return_list$jac %*% w) / c(-1, 2, -3, 4))
    } else {
      wT_gamma <- t(w) %*% X_skew_t_params$gamma
      wT_Sigma_w <- sum((X_skew_t_params$chol_Sigma %*% w) ** 2)

      return_list$jac <- rbind("grad1" = as.vector(X_skew_t_params$mu + X_skew_t_params$a$a11 * X_skew_t_params$gamma),
                               "grad2" = as.vector(2 * X_skew_t_params$a$a21 * X_skew_t_params$scatter %*% w + X_skew_t_params$a$a22 * 2 * as.numeric(wT_gamma) * as.matrix(X_skew_t_params$gamma)),
                               "grad3" = as.vector(X_skew_t_params$a$a31 * 3 * (as.numeric(wT_gamma) ** 2) * as.matrix(X_skew_t_params$gamma) + X_skew_t_params$a$a32 * (as.numeric(t(w) %*% (X_skew_t_params$scatter) %*% w) * as.matrix(X_skew_t_params$gamma) + 2 * as.numeric(wT_gamma) * X_skew_t_params$scatter %*% w )),
                               "grad4" = as.vector(4 * X_skew_t_params$a$a41 * as.numeric((wT_gamma) ** 3) *as.matrix(X_skew_t_params$gamma) + X_skew_t_params$a$a42 * (2 * as.numeric((wT_gamma) ** 2) * X_skew_t_params$scatter %*% w + 2 * as.numeric(wT_Sigma_w) * as.numeric((wT_gamma)) *as.matrix(X_skew_t_params$gamma)) + 4 * X_skew_t_params$a$a43 * (as.numeric(wT_Sigma_w)) * (X_skew_t_params$scatter) %*% w))
      return_list$obj <- sum(lambda * as.vector(return_list$jac %*% w) / c(-1, 2, -3, 4))
    }

    return(return_list)
  }

  # initialization
  start_time <- proc.time()[3]
  if (method == "L-MVSK" || method == "DC") {
    # Define the function to compute several bounds for the second order derivatives
    Compute_bounds <- function(X_skew_t_params) {
      b2 <- 2 * norm(X_skew_t_params$a$a21 * X_skew_t_params$scatter + X_skew_t_params$a$a22 * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma) ,"I")

      H3 <- function(w, X_skew_t_params) {
        H3_mat <-  X_skew_t_params$a$a31 * 6 * as.numeric(t(w) %*% X_skew_t_params$gamma) * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma) + X_skew_t_params$a$a32 * (2 * X_skew_t_params$gamma %*% t(w) %*% X_skew_t_params$scatter + 2 * X_skew_t_params$scatter %*% w %*% t(X_skew_t_params$gamma) + 2 * as.numeric(t(w) %*% X_skew_t_params$gamma) * X_skew_t_params$scatter )
        return(H3_mat)
      }
      H4 <- function(w, X_skew_t_params) {
        H4_mat <- 12 * X_skew_t_params$a$a41 * (as.numeric(t(w) %*% X_skew_t_params$gamma)**2) * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma)  + 2 * X_skew_t_params$a$a42 * (2 * as.numeric(t(w) %*% X_skew_t_params$gamma) * X_skew_t_params$scatter %*% w %*% t(X_skew_t_params$gamma) + ((as.numeric(t(w) %*% X_skew_t_params$gamma))**2) * X_skew_t_params$scatter + 2 * as.numeric(t(w) %*% X_skew_t_params$gamma) * X_skew_t_params$gamma %*% t(w) %*% X_skew_t_params$scatter + as.numeric(t(w) %*% X_skew_t_params$scatter %*% w) * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma)    ) + 4 * X_skew_t_params$a$a43 * (as.numeric(t(w) %*% X_skew_t_params$scatter %*% w) * X_skew_t_params$scatter + 2 * X_skew_t_params$scatter %*% w %*% t(w) %*% X_skew_t_params$scatter)
        return(H4_mat)
      }
      H4_variant <- function(w1, w2, X_skew_t_params) {
        H4_variant_mat <- 12 * X_skew_t_params$a$a41 * (as.numeric(t(w1) %*% X_skew_t_params$gamma) * as.numeric(t(w2) %*% X_skew_t_params$gamma)) * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma)  + 2 * X_skew_t_params$a$a42 * ( as.numeric(t(w1) %*% X_skew_t_params$gamma) * X_skew_t_params$scatter %*% w2 %*% t(X_skew_t_params$gamma) + as.numeric(t(w1) %*% X_skew_t_params$gamma) * X_skew_t_params$scatter %*% w1 %*% t(X_skew_t_params$gamma) + ((as.numeric(t(w1) %*% X_skew_t_params$gamma)) * (as.numeric(t(w2) %*% X_skew_t_params$gamma)) ) * X_skew_t_params$scatter + as.numeric(t(w1) %*% X_skew_t_params$gamma) * X_skew_t_params$gamma %*% t(w2) %*% X_skew_t_params$scatter + as.numeric(t(w2) %*% X_skew_t_params$gamma) * X_skew_t_params$gamma %*% t(w1) %*% X_skew_t_params$scatter + as.numeric(t(w1) %*% X_skew_t_params$scatter %*% w2) * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma)    ) + 4 * X_skew_t_params$a$a43 * (as.numeric(t(w1) %*% X_skew_t_params$scatter %*% w2) * X_skew_t_params$scatter + X_skew_t_params$scatter %*% w1 %*% t(w2) %*% X_skew_t_params$scatter+ X_skew_t_params$scatter %*% w2 %*% t(w1) %*% X_skew_t_params$scatter)
        return(H4_variant_mat)
      }

      w_b <- rep(0, N)
      H3_mat <- matrix(0, N, N)
      for(i in 1:N) {
        w_i <- w_b
        w_i[i] <- 1
        H3_mat <- pmax(H3_mat,  abs(H3(w_i, X_skew_t_params)))
      }
      b3 <- norm(H3_mat, "I")

      H4_mat_variant <- matrix(0, N, N)
      for(i in 1:N) {
        for(j in 1:N) {
          w_i <- w_b
          w_i[i] <- 1
          w_j <- w_b
          w_j[j] <- 1
          H4_mat_variant <- pmax(H4_mat_variant, abs(H4_variant(w_i, w_j, X_skew_t_params)))
        }
      }
      b4 <- norm(H4_mat_variant , "I")
      return(list(b2 = b2, b3 = b3, b4 = b4))
    }

  }
  if (method == "L-MVSK") {
    bounds <- Compute_bounds(X_skew_t_params)
    rho <- lambda[3] * bounds$b3 + lambda[4] * bounds$b4
    H2 <- 2 * (X_skew_t_params$a$a21 * X_skew_t_params$scatter + X_skew_t_params$a$a22 * X_skew_t_params$gamma %*% t(X_skew_t_params$gamma))
  }
  if (method == "DC") {
    bounds <- Compute_bounds(X_skew_t_params)
    rho <- lambda[2]*bounds$b2 + lambda[3]*bounds$b3 + lambda[4]*bounds$b4
  }
  if (method == "PGD" || method == "RFPA" || method == "SQUAREM") {
    # define the PGD update function via water-filling algorithm
    PGD_update <- function(w, eta, g) {
      r <- w - eta * g
      if(sum(r-min(r)) <= 1) {
        gamma_tilde <- (1/N) * (-1+sum(r))
        w <- pmax(0, r - gamma_tilde)
      } else {
        # otherwise we need to find out the position via bisection
        r_vec <- sort(r)
        inv_cum_r_vec <- rev(cumsum(rev(r_vec)))

        n_down <- 0
        n_up <- N
        n <- round((n_down + n_up)/2)

        stop_sign <- FALSE
        while(n_up - n_down > 1) {
          fn <- inv_cum_r_vec[n] - (N-n+1) * r_vec[n]
          if(fn < 1) {
            n_up <- n
          } else if(fn > 1) {
            n_down <- n
          }
          n <- round((n_down + n_up)/2)
        }
        gamma_tilde <- (sum(inv_cum_r_vec[n_up]) -1)/(N-n_up+1)
        w <- pmax(0, r - gamma_tilde)
      }
      return(w)
    }

    # get the gradients from jacobians
    get_gradient <- function(jac) {
      return(-(lambda[1]*jac[1, ] - lambda[2]*jac[2, ] + lambda[3]*jac[3, ] - lambda[4]*jac[4, ]))
    }
  }

  wk <- w_init
  cpu_time <- c(0)
  objs  <- c()
  fun_k <- fun_eval(wk)
  objs <- c(objs, fun_k$obj)


  for (iter in 1:maxiter) {
    # record previous w
    w_old <- wk

    ## construct QP approximation problem (the symbol and scale is adjusted to match the format of solver quadprog::solve.QP)
    switch(method,
           "Q-MVSK" = {
             # PSD approximation
             H_ncvx <- - lambda[3] * fun_k$H3 + lambda[4] * fun_k$H4
             eigen_docom <- eigen(H_ncvx)
             d <- eigen_docom$values
             if(min(d)<0) {
               d[which(d<0)] <- 0
               H_ncvx <- eigen_docom$vectors %*% diag(d) %*% t(eigen_docom$vectors)
             }
             Qk <- H_ncvx + lambda[2] * fun_k$H2  + diag(tau_w, N)

             # solve QP problem
             qk <- lambda[1]*fun_k$jac[1, ] + lambda[3]*fun_k$jac[3, ] - lambda[4]*fun_k$jac[4, ] + H_ncvx%*%wk + tau_w*wk

             sc <- norm(Qk,"2")
             w_hat <- (solve.QP(Qk/sc, qk/sc, t(rbind(rep(1,N), diag(N))),  c(1, rep(0,N)), meq = 1))$solution
             w_hat[which(w_hat<0)] <- 0
             w_hat <- w_hat/sum(w_hat)

             # update w
             wk <- wk + gamma * (w_hat - wk)
             gamma <- gamma * (1 - zeta * gamma)
           },
           "L-MVSK" = {
             Qk <- lambda[2] * H2 + rho * diag(N)
             qk <- rho * wk - (- lambda[1] * (X_skew_t_params$mu + X_skew_t_params$a$a11 * X_skew_t_params$gamma) - lambda[3] * (X_skew_t_params$a$a31 * 3 * (as.numeric(t(wk) %*% X_skew_t_params$gamma) ** 2) * as.matrix(X_skew_t_params$gamma) + X_skew_t_params$a$a32 * (as.numeric(t(wk) %*% (X_skew_t_params$scatter) %*% wk) * as.matrix(X_skew_t_params$gamma) + 2 * as.numeric(t(wk) %*% X_skew_t_params$gamma) * X_skew_t_params$scatter %*% wk )) + lambda[4] * (4 * X_skew_t_params$a$a41 * as.numeric((t(wk) %*% X_skew_t_params$gamma) ** 3) *as.matrix(X_skew_t_params$gamma) + X_skew_t_params$a$a42 * (2 * as.numeric((t(wk) %*% X_skew_t_params$gamma) ** 2) * X_skew_t_params$scatter %*% wk + 2 * as.numeric(t(wk) %*% (X_skew_t_params$scatter) %*% wk) * as.numeric((t(wk) %*% X_skew_t_params$gamma)) *as.matrix(X_skew_t_params$gamma)) + 4 * X_skew_t_params$a$a43 * (as.numeric(t(wk) %*% (X_skew_t_params$scatter) %*% wk)) * (X_skew_t_params$scatter) %*% wk   ))

             wk <- (quadprog::solve.QP(Qk, qk, t(rbind(rep(1,N), diag(N))),  c(1, rep(0,N)), meq = 1))$solution
           },
           "DC" = {
             Qk <- diag(rho, N)
             qk <- rho*wk + lambda[1]*fun_k$jac[1, ] - lambda[2]*fun_k$jac[2, ] + lambda[3]*fun_k$jac[3, ] - lambda[4]*fun_k$jac[4, ]
             sc <- norm(Qk,"2")
             wk <- quadprog::solve.QP(Dmat = Qk/sc, dvec = cbind(qk)/sc, Amat = Amat, bvec = bvec, meq = 1)$solution
           },
           "PGD" ={
             current_obj <- objs[length(objs)]
             # compute the gradient
             eta <- initial_eta
             gk <- get_gradient(fun_k$jac)

             wk_next <- PGD_update(wk, eta, gk)
             fun_k_next <- fun_eval(wk_next)
             next_obj <- fun_k_next$obj

             # backtracking line search
             while(next_obj > current_obj + t(gk) %*% (wk_next - wk) + ((1/(2*eta)) * sum((wk-wk_next)**2)) ){
               eta <- eta * beta
               wk_next <- PGD_update(wk, eta, gk)
               fun_k_next <- fun_eval(wk_next)
               next_obj <- fun_k_next$obj
             }
             wk <- wk_next
           },
           "RFPA" = {
             current_obj <- objs[length(objs)]

             # Try SQUAREM acceleration
             gk <- get_gradient(fun_k$jac)

             # One iterate of update
             wk1 <- PGD_update(wk, 1/tau, gk)
             fun_k1 <- fun_eval(wk1)

             # Another iterate of update
             wk2 <- PGD_update(wk1, 1/tau, get_gradient(fun_k1$jac) )

             r <- wk1 - wk
             v <- (wk2 - wk1) - r
             if(max(abs(v)) == 0) {
               alpha <- -1
             } else {
               alpha <- - (norm(r, "2"))/(norm(v, "2"))
             }

             if(as.numeric(t(r) %*% v) < 0)
             {
               alpha <- max(alpha, (as.numeric(t(r) %*% r))/(as.numeric(t(r) %*% v)))
             }

             wkt <- wk - 2 * alpha * r + alpha * alpha * v
             fun_kt <- fun_eval(wkt)

             wk_next <- PGD_update(wkt, 1/tau, get_gradient(fun_kt$jac) )
             #wk_next <- PGD_update(wkt, 1/tau, get_gradient(fun_kt$jac) )
             fun_k_next <- fun_eval(wk_next)

             next_obj <- fun_k_next$obj

             # If we need PGD to leave the current point
             if(next_obj > current_obj) {
               eta <- initial_eta

               wk_next <- PGD_update(wk, eta, gk)
               fun_k_next <- fun_eval(wk_next)
               next_obj <- fun_k_next$obj

               # backtracking line search
               while(next_obj > current_obj + t(gk) %*% (wk_next - wk) + ((1/(2*eta)) * sum((wk-wk_next)**2)) ){
                 eta <- eta * beta
                 wk_next <- PGD_update(wk, eta, gk)
                 fun_k_next <- fun_eval(wk_next)
                 next_obj <- fun_k_next$obj
               }
               #print(eta)
             }
             wk <- wk_next
           },
           "SQUAREM" = {
             current_obj <- objs[length(objs)]

             # Try SQUAREM acceleration
             gk <- get_gradient(fun_k$jac)

             # One iterate of update
             wk1 <- PGD_update(wk, 1/tau, gk)
             fun_k1 <- fun_eval(wk1)

             # Another iterate of update
             wk2 <- PGD_update(wk1, 1/tau, get_gradient(fun_k1$jac) )

             r <- wk1 - wk
             v <- (wk2 - wk1) - r
             if(max(abs(v)) == 0) {
               alpha <- -1
             } else {
               alpha <- min(-1, - (norm(r, "2"))/(norm(v, "2")))
             }

             wkt <- wk - 2 * alpha * r + alpha * alpha * v
             fun_kt <- fun_eval(wkt)

             wk_next <- PGD_update(wkt, 1/tau, get_gradient(fun_kt$jac) )
             fun_k_next <- fun_eval(wk_next)

             next_obj <- fun_k_next$obj

             # backtracking line search on alpha
             while(next_obj > current_obj) {
               alpha <- 0.5 * (alpha - 1)
               wkt <- wk - 2 * alpha * r + alpha * alpha * v
               fun_kt <- fun_eval(wkt)

               wk_next <- PGD_update(wkt, 1/tau, get_gradient(fun_kt$jac) )
               fun_k_next <- fun_eval(wk_next)

               next_obj <- fun_k_next$obj
             }
             wk <- wk_next
           },
           stop("Method unknown")
    )

    wk[which(wk<0)] <- 0
    wk <- wk/sum(wk)

    # recording...
    cpu_time <- c(cpu_time, proc.time()[3] - start_time)

    if(method == "PGD" || method == "RFPA") {
      fun_k <- fun_k_next
    } else {
      fun_k <- fun_eval(wk)
    }

    objs <- c(objs, fun_k$obj)

    # termination criterion
    # has_w_converged <- all(abs(wk - w_old) <= .5 * wtol )
    # has_w_converged <- norm(wk - w_old, "2") <= wtol * norm(w_old, "2")
    has_w_converged <- sqrt(sum((wk - w_old)**2))/ sqrt(sum((w_old)**2)) < wtol
    has_f_converged <- abs(diff(tail(objs, 2))) < ftol
    has_cross_stopval <- tail(objs, 1) <= stopval

    if(is.infinite(stopval)) {
      if (has_w_converged && has_f_converged) break
    } else {
      if(has_cross_stopval) break
    }
  }

  return(list(
    "w"                      = wk,
    "cpu_time_vs_iterations" = cpu_time,
    "objfun_vs_iterations"   = objs * lambda_max,
    "iterations"             = 0:iter,
    "convergence"            = !(iter == maxiter),
    "moments"                = as.vector(fun_k$jac %*% wk) / c(1, 2, 3, 4)
  ))
}
