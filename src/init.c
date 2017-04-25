#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GDINA_aggregateCol(SEXP, SEXP);
extern SEXP GDINA_designM(SEXP, SEXP);
extern SEXP GDINA_fitstats(SEXP, SEXP, SEXP);
extern SEXP GDINA_Lik(SEXP, SEXP, SEXP, SEXP);
extern SEXP GDINA_NgRg(SEXP, SEXP, SEXP);
extern SEXP GDINA_Rljs(SEXP, SEXP, SEXP);
extern SEXP GDINA_scorefun(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GDINA_sequP(SEXP, SEXP, SEXP);
extern SEXP GDINA_uP(SEXP, SEXP);
extern SEXP GDINA_varsigma(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GDINA_aggregateCol", (DL_FUNC) &GDINA_aggregateCol, 2},
    {"GDINA_designM",      (DL_FUNC) &GDINA_designM,      2},
    {"GDINA_fitstats",     (DL_FUNC) &GDINA_fitstats,     3},
    {"GDINA_Lik",          (DL_FUNC) &GDINA_Lik,          4},
    {"GDINA_NgRg",         (DL_FUNC) &GDINA_NgRg,         3},
    {"GDINA_Rljs",         (DL_FUNC) &GDINA_Rljs,         3},
    {"GDINA_scorefun",     (DL_FUNC) &GDINA_scorefun,     5},
    {"GDINA_sequP",        (DL_FUNC) &GDINA_sequP,        3},
    {"GDINA_uP",           (DL_FUNC) &GDINA_uP,           2},
    {"GDINA_varsigma",     (DL_FUNC) &GDINA_varsigma,     3},
    {NULL, NULL, 0}
};

void R_init_GDINA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
