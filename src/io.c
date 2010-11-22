#include <R.h>
#include <Rinternals.h>

SEXP convert_quality(SEXP ascii_quality) {

  if (!isString(ascii_quality))
    error("ascii_quality must be of type 'character'");

  R_len_t l, i;
  l = length(ascii_quality);
  for (i=0; i<=l; i++) {
    printf("%d", CHAR(STRING_ELT(ascii_quality, i)));
  }
  return NILSXP;
}
