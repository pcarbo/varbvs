#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "varbvsr.h"

// See "Registering native routines" in "Writing R Extensions" manual
// for an explanation of what these lines of code do.

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

const static R_CallMethodDef R_CallDef[] = {
  CALLDEF(diagsq_Call,3),
  CALLDEF(diagsqt_Call,3),
  CALLDEF(varbvsnormupdate_Call,10),
  CALLDEF(varbvsbinupdate_Call,11),
  CALLDEF(varbvsbinzupdate_Call,11),
  CALLDEF(varbvs_varbvsnormupdate_rcpp,10),
  {NULL, NULL, 0}
};

void attribute_visible R_init_varbvs(DllInfo *dll)
{
  R_registerRoutines(dll,NULL,R_CallDef,NULL,NULL);
  R_useDynamicSymbols(dll,FALSE);
  R_forceSymbols(dll,TRUE);
}
