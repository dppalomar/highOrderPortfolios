// This file is a copy of several internal functions in R package PerformanceAnalytics (version: 2.0.4)
// for full details see: https://cran.r-project.org/web/packages/PerformanceAnalytics/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP M3mat2vec(SEXP, SEXP);
extern SEXP M3vec2mat(SEXP, SEXP);
extern SEXP M4mat2vec(SEXP, SEXP);
extern SEXP M4vec2mat(SEXP, SEXP);
extern SEXP M3port_grad(SEXP, SEXP, SEXP);
extern SEXP M3port(SEXP, SEXP, SEXP);
extern SEXP M4port(SEXP, SEXP, SEXP);
extern SEXP M4port_grad(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"M3mat2vec",        (DL_FUNC) &M3mat2vec,         2},
  {"M3vec2mat",        (DL_FUNC) &M3vec2mat,         2},
  {"M4mat2vec",        (DL_FUNC) &M4mat2vec,         2},
  {"M4vec2mat",        (DL_FUNC) &M4vec2mat,         2},
  {"M3port_grad",      (DL_FUNC) &M3port_grad,       3},
  {"M3port",           (DL_FUNC) &M3port,            3},
  {"M4port",           (DL_FUNC) &M4port,            3},
  {"M4port_grad",      (DL_FUNC) &M4port_grad,       3},
  {NULL, NULL, 0}
};

void R_init_highOrderPortfolios(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}