#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>

#define LINE_BUFFER 500

static const int LINES_PER_FASTQ_REC = 4;

typedef struct {
  char header[LINE_BUFFER];
  char sequence[LINE_BUFFER];
  char quality[LINE_BUFFER];
} fastq_block;

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

int getline(FILE *fp, char *buf, int bufsize) {
  /* 
     Ported from Rconnections.h, because this hasn't been opened up to
     package developers yet. This verison works with a raw C stream
     rather than an R connection.
  */
  int c, nbuf = -1;
  
  while((c = fgetc(fp)) != EOF) {
    if (nbuf+1 >= bufsize) error("Line longer than buffer size");
    if (c != '\n') {
      buf[++nbuf] = c;
    } else {
      buf[++nbuf] = '\0';
      break;
    }
  }
  /* Make sure it is null-terminated and count is correct, even if
   *  file did not end with newline.
   */
  if (nbuf >= 0 && buf[nbuf]) {
    if (nbuf+1 >= bufsize) error("Line longer than buffer size");
    buf[++nbuf] = '\0';
  }
  return nbuf;
}

char *trim(char *str) {
  char *end;
  while (isspace(*str)) str++;

  if (*str == '\0')
    return str;

  end = str + strlen(str) - 1;
  while (end > str && isspace(*end)) end--;

  *(end+1) = 0;
  return str;
}

fastq_block *_get_fastq_block(FILE *fp) {
  unsigned int nlines = 0, i;
  fastq_block *block = malloc(sizeof(fastq_block));

  if (fp == NULL) error("stream pointer is null");
  
  char *buf = malloc(LINE_BUFFER);
  if (!buf) error("cannot allocate buffer in read_fastq_file");
  
  for (i = 0; i < 4; i++) {
    if (getline(fp, buf, LINE_BUFFER) == -1) return NULL;

    switch (nlines % LINES_PER_FASTQ_REC) {
    case 0:
    case 2:
      buf++; // ditch header character
      strcpy(block->header, trim(buf));
      break;
    case 1:
      // sequence
      strcpy(block->sequence, trim(buf));
      break;
    case 3:
      // quality
      strcpy(block->quality, trim(buf));
      break;
    default:
      error("unexpected error; consult maintainer - %d", nlines);
    }
    
    nlines++;
  }
  free(buf);
  return block;
}

/* SEXP update_sequence_counts(char *seq, SEXP base_counts) { */
/*   /\* */
/*     Given a sequence, update the counts in `base_counts`. */
/*   *\/ */
  
/*   unsigned int i, nseq = strlen(seq); */
/*   SEXP tmp; */
/*   PROTECT(tmp = allocVector(INTSXP, 5)); */
  
/*   for (i = 0; i < nseq; i++) { */
/*     // TODO check length of base_counts */
/*     switch (seq[i]) { */
/*     case 'A': */
/*       tmp = INTEGER(VECTOR_ELT(base_counts, i)); */
/*       printf("-->%d\n", tmp[0]); */
/*       tmp[0]++; */
/*       SET_VECTOR_ELT(base_counts, i, tmp); */
/*       break; */
/*     default: */
/*       break; */
/*     } */
/*   } */
/*   UNPROTECT(1); */
/*   return base_counts; */
/* } */

SEXP summarize_fastq_file(SEXP filename) {
  long unsigned int nblock = 0;
  fastq_block *tmp;
  SEXP base_counts;
  int *ibc_matrix, i, j;
  
  FILE *fp = fopen(CHAR(STRING_ELT(filename, 0)), "r");
  if (fp == NULL)
    error("failed to open file '%s'", CHAR(STRING_ELT(filename, 0)));

  PROTECT(base_counts = allocMatrix(INTSXP, 5, 100));
  ibc_matrix = INTEGER(base_counts);
  
  int nx = 5, ny = 100;
  for (i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++)
      ibc_matrix[i + nx*j] = 0;
  }

  while ((tmp = _get_fastq_block(fp)) != NULL) {
    for (i = 0; i < strlen(tmp->sequence); i++) {
      // TODO - we could cast these to char, and use a lookup table
      if (tmp->sequence[i] == 'A')
        ibc_matrix[5*i]++;
      if (tmp->sequence[i] == 'C')
        ibc_matrix[5*i + 1]++;
      if (tmp->sequence[i] == 'T')
        ibc_matrix[5*i + 2]++;
      if (tmp->sequence[i] == 'G')
        ibc_matrix[5*i + 3]++;
      if (tmp->sequence[i] == 'N')
        ibc_matrix[5*i + 4]++;
    }

    /* printf("head: %s\n", tmp->header); */
    //printf("seq: %s\n", tmp->sequence);
    /* printf("qual: %s\n", tmp->quality);     */
  }

  UNPROTECT(1);
  return base_counts;
}
