#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "khash.h"

KHASH_MAP_INIT_STR(str, int)

#define LINE_BUFFER 500

#define NUM_BASES 5 // includes N

typedef enum {
  SANGER,
  SOLEXA,
  ILLUMINA
} quality_type;

#define Q_OFFSET 0
#define Q_MIN 1
#define Q_MAX 2

static const int quality_contants[3][3] = {
  // offset, min, max
  {33, 0, 93}, // SANGER
  {64, -5, 62}, // SOLEXA
  {64, 0, 62}, // ILLUMINA
};

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

int get_line(FILE *fp, char *buf, int bufsize) {
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
    if (get_line(fp, buf, LINE_BUFFER) == -1) return NULL;

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

void update_summary_matrices(fastq_block *block, int *base_matrix, int *qual_matrix, quality_type q_type) {
  /*
    Given `fastq_block`, adjust the nucloeotide frequency
    `counts_matrix` accordingly.
  */
  int i, len;
  int q_range = quality_contants[q_type][Q_MAX] - quality_contants[q_type][Q_MIN];
  int q_min = quality_contants[q_type][Q_MIN];
  int q_max = quality_contants[q_type][Q_MAX];
  int q_offset = quality_contants[q_type][Q_OFFSET];
  
  if (strlen(block->sequence) != strlen(block->quality))
    error("improperly formatted FASTQ file; sequence and quality lengths differ");

  len = strlen(block->sequence);

  for (i = 0; i < len; i++) {
    switch ((char) block->sequence[i]) {
    case 'A':
      base_matrix[5*i]++;
      break;
    case 'T':
      base_matrix[5*i + 1]++;
      break;
    case 'C':
      base_matrix[5*i + 2]++;
      break;
    case 'G':
      base_matrix[5*i + 3]++;
      break;
    case 'N':
      base_matrix[5*i + 4]++;
      break;
    default:
      error("Sequence character encountered that is not"
            " A, T, C, G, or N: '%c'", block->sequence[i]);
    }
    
    if ((char) block->quality[i] - q_offset < q_min || (char) block->quality[i] - q_offset > q_max)
      error("base quality out of range (%d < b < %d) encountered: %d", q_min, 
            q_max, (char) block->quality[i]);
    
    qual_matrix[(q_range+1)*i + ((char) block->quality[i]) - q_offset - q_min]++;
  }
}

void zero_int_matrix(int *matrix, int nx, int ny) {
  int i, j;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++)
      matrix[i + nx*j] = 0;
  }
}

SEXP summarize_fastq_file(SEXP filename, SEXP max_length, SEXP quality_type, SEXP hash, SEXP verbose) {
  if (!isString(filename))
    error("filename should be a string");
  if (INTEGER(max_length)[0] > LINE_BUFFER)
    error("You have specified a max_length less than the C buffer size. "
          "Adjust the 'LINE_BUFFER' macro and recompile to run sequences"
          " greater than %d.", LINE_BUFFER);

  khash_t(str) *h;
  khiter_t k;
  int is_missing, ret, size_out_list = 3;
  unsigned int num_unique_seqs = 0;

  unsigned int nblock = 0;
  fastq_block *block;
  SEXP base_counts, qual_counts, seq_hash, seq_hash_names, out_list, seq_lengths;
  int *ibc, *iqc, *isl, i, j, q_type, q_range;

  q_type = INTEGER(quality_type)[0];
  q_range = quality_contants[q_type][Q_MAX] - quality_contants[q_type][Q_MIN];
  
  if (LOGICAL(hash)[0]) {
    h = kh_init(str);
    kh_resize(str, h, 1572869);
    size_out_list = 4;
  }

  FILE *fp = fopen(CHAR(STRING_ELT(filename, 0)), "r");
  if (fp == NULL)
    error("failed to open file '%s'", CHAR(STRING_ELT(filename, 0)));

  PROTECT(out_list = allocVector(VECSXP, size_out_list));
  PROTECT(base_counts = allocMatrix(INTSXP, NUM_BASES, INTEGER(max_length)[0]));
  PROTECT(qual_counts = allocMatrix(INTSXP, q_range + 1, INTEGER(max_length)[0]));
  PROTECT(seq_lengths = allocVector(INTSXP, INTEGER(max_length)[0]));
  
  ibc = INTEGER(base_counts);
  iqc = INTEGER(qual_counts);
  isl = INTEGER(seq_lengths);

  zero_int_matrix(ibc, NUM_BASES, INTEGER(max_length)[0]);
  zero_int_matrix(iqc, q_range + 1, INTEGER(max_length)[0]);

  while ((block = read_fastq_block(fp)) != NULL) {
    R_CheckUserInterrupt();
     
    update_summary_matrices(block, ibc, iqc, q_type);
    
    isl[nblock] = strlen(block->sequence);

    if (LOGICAL(hash)[0]) {
      k = kh_get(str, h, block->sequence);
      is_missing = (k == kh_end(h));
      if (is_missing) {
        k = kh_put(str, h, strdup(block->sequence), &ret);
        kh_value(h, k) = 1;
        num_unique_seqs++;
      } else 
        kh_value(h, k) = kh_value(h, k) + 1;
      
      if (LOGICAL(verbose)[0] && nblock % 100000 == 0)
        printf("on block %d, %d entries in hash table\n", nblock, num_unique_seqs);
    }
    free(block);
    nblock++;
  }

  PROTECT(seq_hash = allocVector(VECSXP, num_unique_seqs));
  PROTECT(seq_hash_names = allocVector(VECSXP, num_unique_seqs));

  if (LOGICAL(verbose)[0])
    printf("processing complete... now loading C hash structure to R...\n");

  i = 0;
  for (k = kh_begin(h); k != kh_end(h); ++k) {
    if (kh_exist(h, k)) {
      SET_VECTOR_ELT(seq_hash_names, i, mkString(kh_key(h, k)));
      SET_VECTOR_ELT(seq_hash, i, ScalarInteger(kh_value(h, k)));
      /* per the comment here
         (http://attractivechaos.wordpress.com/2009/09/29/khash-h/),
         using character arrays keys with strdup must be freed during
         table traverse. */
      free((char *) kh_key(h, k));
      i++;
    }
  }
    
  kh_destroy(str, h);

  SET_VECTOR_ELT(out_list, 0, base_counts);
  SET_VECTOR_ELT(out_list, 1, qual_counts);
  SET_VECTOR_ELT(out_list, 2, seq_lengths);
  if (LOGICAL(hash)[0]) {
    setAttrib(seq_hash, R_NamesSymbol, seq_hash_names);
    SET_VECTOR_ELT(out_list, 3, seq_hash);
  }

  UNPROTECT(6);
  fclose(fp);
  return out_list;
}
