##
## User installation
##
# Local installation
install.packages(file.choose(), repos = NULL, type = "source")
# Installation from GitHub
devtools::install_github("dppalomar/highOrderPortfolios")
# Installation from CRAN
install.packages("highOrderPortfolios")
# Getting help
library(highOrderPortfolios)
help(package = "highOrderPortfolios")
package?highOrderPortfolios
?design_MVSK_portfolio
citation("highOrderPortfolios")
vignette(package = "highOrderPortfolios")


##
## Developer commands (https://r-pkgs.org/)
##
devtools::load_all()  #or Ctrl-Shift-L
devtools::document()  #to generate all documentation via roxygen
devtools::install()
devtools::install(dependencies = FALSE)
library(highOrderPortfolios)
#tools::showNonASCIIfile("R/MVSK.R")


# Code tests
devtools::test()
#covr::package_coverage()  #coverage of tests


# CRAN check and submission (https://r-pkgs.org/release.html)
#  checklist: https://kalimu.github.io/post/checklist-for-r-package-submission-to-cran/
devtools::check()  # run_dont_test = TRUE
rcmdcheck::rcmdcheck()  # build_args = "--run-donttest"
devtools::build()
#devtools::revdep(pkg = "highOrderPortfolios")  # to check reverse dependencies
#devtools::check_win_release()  #to check under windows
#R CMD build . --resave-data  # this is to generate tarball
#R CMD check highOrderPortfolios_0.0.2.tar.gz --as-cran --run-donttest  # this is before submission to CRAN
#R CMD install highOrderPortfolios_0.0.2.tar.gz
#submit the tarball directly via the webform: https://cran.r-project.org/submit.html
