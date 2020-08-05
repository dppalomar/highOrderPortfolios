# This file is a copy of several internal functions in R package PerformanceAnalytics (version: 2.0.4)
# for full details see: https://cran.r-project.org/web/packages/PerformanceAnalytics/
# Copyright (c) 2004-2020 Kris Boudt and Brian G. Peterson
# Copyright (c) 2004-2020 Peter Carl and Brian G. Peterson for PerformanceAnalytics


# Wrapper function for casting the vector with unique coskewness elements into the matrix format
PerformanceAnalytics_M3.vec2mat <- function(M3, p) {
  # M3        : numeric vector with unique coskewness elements (p * (p + 1) * (p + 2) / 6)
  # p         : dimension of the data
  
  if (NCOL(M3) != 1) stop("M3 must be a vector")
  
  .Call('M3vec2mat', M3, as.integer(p), PACKAGE = "highOrderPortfolios")
}

# Wrapper function for casting the coskewness matrix into the vector of unique elements
PerformanceAnalytics_M3.mat2vec <- function(M3) {
  # M3        : numeric matrix of dimension p x p^2
  
  if (is.null(dim(M3))) stop("M3 must be a matrix")
  
  .Call('M3mat2vec', as.numeric(M3), NROW(M3), PACKAGE = "highOrderPortfolios")
}

# Wrapper function for casting the cokurtosis matrix into the vector of unique elements
PerformanceAnalytics_M4.mat2vec <- function(M4) {
  # M4        : numeric matrix of dimension p x p^3
  
  if (is.null(dim(M4))) stop("M4 must be a matrix")
  
  .Call('M4mat2vec', as.numeric(M4), NROW(M4), PACKAGE="highOrderPortfolios")
}

# Wrapper function for casting the vector with unique cokurtosis elements into the matrix format
PerformanceAnalytics_M4.vec2mat <- function(M4, p) {
  # M4        : numeric vector with unique cokurtosis elements (p * (p + 1) * (p + 2) * (p + 3) / 24)
  # p         : dimension of the data
  
  if (NCOL(M4) != 1) stop("M4 must be a vector")
  
  .Call('M4vec2mat', M4, as.integer(p), PACKAGE = "highOrderPortfolios")
}


PerformanceAnalytics_portm3 <- function(w, M3)
{
  w <- as.numeric(w)
  if (NCOL(M3) != 1) M3 <- PerformanceAnalytics_M3.mat2vec(M3)  # M3 <- M3.mat2vec(M3)
  .Call('M3port', w, M3, length(w), PACKAGE = "highOrderPortfolios")
}

PerformanceAnalytics_derportm3 <- function(w, M3)
{
  w <- as.numeric(w)
  if (NCOL(M3) != 1) M3 <- PerformanceAnalytics_M3.mat2vec(M3) #  M3 <- M3.mat2vec(M3)
  as.matrix(.Call('M3port_grad', w, M3, length(w), PACKAGE = "highOrderPortfolios"), ncol = 1)
}

PerformanceAnalytics_portm4 <- function(w, M4)
{
  w <- as.numeric(w)
  if (NCOL(M4) != 1) M4 <- PerformanceAnalytics_M4.mat2vec(M4)  # M4 <- M4.mat2vec(M4)
  .Call('M4port', w, M4, length(w), PACKAGE = "highOrderPortfolios")
}

PerformanceAnalytics_derportm4 <- function(w, M4)
{
  w <- as.numeric(w)
  if (NCOL(M4) != 1) M4 <- PerformanceAnalytics_M4.mat2vec(M4)  # M4 <- M4.mat2vec(M4)
  as.matrix(.Call('M4port_grad', w, M4, length(w), PACKAGE = "highOrderPortfolios"), ncol = 1)
}
