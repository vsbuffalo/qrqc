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

fastq_block *read_fastq_block(FILE *fp) {
  unsigned int nlines = 0, i;
  fastq_block *block = malloc(sizeof(fastq_block));

  if (fp == NULL) error("stream pointer is null");
  
  char *obuf, *buf = malloc(LINE_BUFFER);
  obuf = buf;
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
  free(obuf);
  return block;
}

void update_summary_matrices(fastq_block *block, int *base_matrix, int *qual_matrix) {
  /*
    Given `fastq_block`, adjust the nucloeotide frequency
    `counts_matrix` accordingly.
  */
  int i, len;
  
  if (strlen(block->sequence) != strlen(block->quality))
    error("improperly formatted FASTQ file; sequence and quality lengths differ");

  len = strlen(block->sequence);

  for (i = 0; i < len; i++) {
    switch ((int) block->sequence[i]) {
    case 'A':
      base_matrix[5*i]++;
      break;
    case 'C':
      base_matrix[5*i + 1]++;
      break;
    case 'T':
      base_matrix[5*i + 2]++;
      break;
    case 'G':
      base_matrix[5*i + 3]++;
      break;
    case 'N':
      base_matrix[5*i + 4]++;
    default:
      error("Sequence character encountered that is not A, T, C, G, or N");
    }
  }
}

void zero_matrix(int *matrix, int nx, int ny) {
  int i, j;
  for (i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++)
      matrix[i + nx*j] = 0;
  }
}

SEXP summarize_fastq_file(SEXP filename) {
  long unsigned int nblock = 0;
  fastq_block *block;
  SEXP base_counts, qual_counts;
  int *ibc, *iqc, i, j, nx = 5, ny = 100; // TODO make macros/check & warn
  
  FILE *fp = fopen(CHAR(STRING_ELT(filename, 0)), "r");
  if (fp == NULL)
    error("failed to open file '%s'", CHAR(STRING_ELT(filename, 0)));

  PROTECT(base_counts = allocMatrix(INTSXP, 5, 100));
  PROTECT(qual_counts = allocMatrix(INTSXP, 5, 100));
  
  ibc = INTEGER(base_counts);
  iqc = INTEGER(qual_counts);
  
  zero_matrix(ibc, nx, ny);

  while ((block = read_fastq_block(fp)) != NULL) {
    update_summary_matrices(block, ibc, iqc);
    free(block);
  }

  UNPROTECT(2);
  return base_counts;
}
