#include <R.h>
#include <Rinternals.h>
#include <string.h>

SEXP string_to_base_qualities(SEXP ascii_quality) {

  if (!isString(ascii_quality))
    error("ascii_quality must be of type 'character'");  

  SEXP quals;
  int i, l = strlen(CHAR(STRING_ELT(ascii_quality, 0)));

  PROTECT(quals = allocVector(INTSXP, l));

  for (i = 0; i < l; i++) {
    //printf("base: %c\n", CHAR(STRING_ELT(ascii_quality, 0))[i]);
    INTEGER(quals)[i] = (int) CHAR(STRING_ELT(ascii_quality, 0))[i] - 33;
  }
  UNPROTECT(1);
  return quals;
}


