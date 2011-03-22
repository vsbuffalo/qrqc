#include "io.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
  {"summarize_file", (DL_FUNC) &summarize_file, 5},
  {NULL, NULL, 0}
};

void R_init_io(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
