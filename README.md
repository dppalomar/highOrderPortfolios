
<!-- README.md is generated from README.Rmd. Please edit that file -->





# highOrderPortfolios
<!---
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/highOrderPortfolios)](https://CRAN.R-project.org/package=highOrderPortfolios)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/highOrderPortfolios)](https://CRAN.R-project.org/package=highOrderPortfolios)
[![CRAN Downloads Total](https://cranlogs.r-pkg.org/badges/grand-total/highOrderPortfolios?color=brightgreen)](https://CRAN.R-project.org/package=highOrderPortfolios)
--->

The classical Markowitz's mean-variance portfolio formulation ignores 
heavy tails and skewness. High-order portfolios use higher order moments to
better characterize the return distribution. Different formulations and fast 
algorithms are proposed for high-order portfolios based on the mean, variance, 
skewness, and kurtosis.
The package is based on the papers [Zhou and Palomar (2021)](https://arxiv.org/abs/2008.00863)
and [Wang, Zhou, Ying, and Palomar (2022)](https://arxiv.org/pdf/2206.02412v1.pdf).


## Installation
The package can be installed from [CRAN](https://CRAN.R-project.org/package=highOrderPortfolios) or [GitHub](https://github.com/dppalomar/highOrderPortfolios):

```r
# install stable version from CRAN
install.packages("highOrderPortfolios")

# install development version from GitHub
devtools::install_github("dppalomar/highOrderPortfolios")
```

To get help:

```r
library(highOrderPortfolios)
help(package = "highOrderPortfolios")
?design_MVSK_portfolio_via_sample_moments
?design_MVSK_portfolio_via_skew_t
```

To cite `highOrderPortfolios` in publications:

```r
citation("highOrderPortfolios")
```


## Usage

```r
library(highOrderPortfolios)
data(X50)

# non-parametric case: estimate sample moments
X_moments <- estimate_sample_moments(X50)

# parametric case: estimate the multivariate skew t distribution
X_skew_t_params <- estimate_skew_t(X50)

# choose hyper-parameter moment weights for the MVSK formulation
xi <- 10
lmd <- c(1, xi/2, xi*(xi+1)/6, xi*(xi+1)*(xi+2)/24)

# design portfolio
sol_nonparam <- design_MVSK_portfolio_via_sample_moments(lmd, X_moments)
sol_param    <- design_MVSK_portfolio_via_skew_t(lmd, X_skew_t_params)

# plot
barplot(cbind("via sample moments" = sol_nonparam$w, "via skew t modeling" = sol_param$w), beside = TRUE, 
        col = c(rep("darkblue", 50), rep("darkred", 50)),
        main = "MSVK portfolio", xlab = "asset indexes", ylab = "portfolio weights")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="75%" style="display: block; margin: auto;" />



## Documentation
For more detailed information, please check the
[vignette](https://htmlpreview.github.io/?https://github.com/dppalomar/highOrderPortfolios/blob/master/vignettes/DesignOfHighOrderPortfolios.html) and  [vignette](https://CRAN.R-project.org/package=highOrderPortfolios/vignettes/highOrderPortfolios.html).



## Links
Package: [CRAN](https://CRAN.R-project.org/package=highOrderPortfolios) and [GitHub](https://github.com/dppalomar/highOrderPortfolios).

README file: [GitHub-readme](https://github.com/dppalomar/highOrderPortfolios/blob/master/README.md).

Vignette: [GitHub-vignette](https://htmlpreview.github.io/?https://github.com/dppalomar/highOrderPortfolios/blob/master/vignettes/DesignOfHighOrderPortfolios.html) and [CRAN-vignette](https://CRAN.R-project.org/package=highOrderPortfolios/vignettes/highOrderPortfolios.html).

