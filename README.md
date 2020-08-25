
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
The package is based on the paper [Zhou and Palomar (2020)](https://arxiv.org/abs/2008.00863).


## Installation
The package can be installed from [GitHub](https://github.com/dppalomar/highOrderPortfolios):
<!---[CRAN](https://CRAN.R-project.org/package=highOrderPortfolios) or---> 

```r
# install stable version from CRAN (not available yet)
#install.packages("highOrderPortfolios")

# install development version from GitHub
devtools::install_github("dppalomar/highOrderPortfolios")
```

To get help:

```r
library(highOrderPortfolios)
help(package = "highOrderPortfolios")
?design_MVSK_portfolio
```

To cite `highOrderPortfolios` in publications:

```r
citation("highOrderPortfolios")
```


## Usage

```r
library(highOrderPortfolios)
data(X50)

# estimate parameters
X_moments <- estimate_moments(X50, adjust_magnitude = TRUE)

# decide problem setting
w0 <- rep(1/50, 50)
w0_moments <- eval_portfolio_moments(w0, X_moments)
d <- abs(w0_moments) 
kappa <- 0.3 * sqrt(w0%*%X_moments$Sgm%*%w0)

# portfolio optimization
sol <- design_MVSKtilting_portfolio(d, X_moments, w_init = w0, w0 = w0, w0_moments = w0_moments, kappa = kappa)

# plot
plot(sol$cpu_time, sol$objs, "b")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="75%" style="display: block; margin: auto;" />

```r
barplot(sol$w)
```

<img src="man/figures/README-unnamed-chunk-5-2.png" width="75%" style="display: block; margin: auto;" />



## Documentation
For more detailed information, please check the
[vignette](https://htmlpreview.github.io/?https://github.com/dppalomar/highOrderPortfolios/blob/master/vignettes/DesignOfHighOrderPortfolios.html).

<!---
[vignette](https://CRAN.R-project.org/package=highOrderPortfolios/vignettes/highOrderPortfolios.html).
--->


## Links
Package: [GitHub](https://github.com/dppalomar/highOrderPortfolios).
<!---[CRAN](https://CRAN.R-project.org/package=highOrderPortfolios) and --->

README file: [GitHub-readme](https://github.com/dppalomar/highOrderPortfolios/blob/master/README.md).

Vignette: [GitHub-vignette](https://htmlpreview.github.io/?https://github.com/dppalomar/highOrderPortfolios/blob/master/vignettes/DesignOfHighOrderPortfolios.html).
<!---[CRAN-vignette](https://CRAN.R-project.org/package=highOrderPortfolios/vignettes/highOrderPortfolios.html).--->

