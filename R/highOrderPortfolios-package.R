#' highOrderPortfolios: Design of High-Order Portfolios via Mean, Variance, Skewness, and Kurtosis
#' 
#' The classical Markowitz's mean-variance portfolio formulation ignores 
#' heavy tails and skewness. High-order portfolios use higher order moments to
#' better characterize the return distribution. Different formulations and fast 
#' algorithms are proposed for high-order portfolios based on the mean, variance, 
#' skewness, and kurtosis.
#'
#'
#' @section Functions:
#' \code{\link{design_MVSK_portfolio_via_sample_moments}()}, \code{\link{design_MVSK_portfolio_via_skew_t}()},
#' and \code{\link{design_MVSKtilting_portfolio_via_sample_moments}()}
#'
#' @section Help:
#' For a quick help see the README file:
#' \href{https://github.com/dppalomar/highOrderPortfolios/blob/master/README.md}{GitHub-README}.
#'
# For more details see the vignette:
# \href{https://CRAN.R-project.org/package=highOrderPortfolios/vignettes/highOrderPortfolios.html}{CRAN-vignette}.
#'
#' @author Rui Zhou, Xiwen Wang, and Daniel P. Palomar
#'
#' @references
#' R. Zhou and D. P. Palomar, "Solving High-Order Portfolios via Successive Convex Approximation Algorithms," 
#' in \emph{IEEE Transactions on Signal Processing}, vol. 69, pp. 892-904, 2021.
#' <https://doi.org/10.1109/TSP.2021.3051369>.
#' 
#' X. Wang, R. Zhou, J. Ying, and D. P. Palomar, "Efficient and Scalable High-Order Portfolios Design via Parametric Skew-t Distribution," 
#' Available in arXiv, 2022. <https://arxiv.org/pdf/2206.02412.pdf>.
#' 
#'
#' @useDynLib highOrderPortfolios
#' @docType package
#' @name highOrderPortfolios-package
NULL
