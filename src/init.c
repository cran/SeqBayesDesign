#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _SeqBayesDesign_TMCMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SeqBayesDesign_TMCMCALT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_SeqBayesDesign_TMCMC",    (DL_FUNC) &_SeqBayesDesign_TMCMC,    9},
  {"_SeqBayesDesign_TMCMCALT", (DL_FUNC) &_SeqBayesDesign_TMCMCALT, 9},
  {NULL, NULL, 0}
};

void R_init_SeqBayesDesign(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
