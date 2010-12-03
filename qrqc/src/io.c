#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "khash.h"
#include "kseq.h"

KHASH_MAP_INIT_STR(str, int)

KSEQ_INIT(gzFile, gzread)

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
  {64, 0, 62} // ILLUMINA
};


void update_base_matrices(kseq_t *block, int *base_matrix) {
  /*
    Given `fastx_block`, adjust the nucloeotide frequency
    `counts_matrix` accordingly.
  */
  int i, len;
  
  len = strlen(block->seq.s);

  for (i = 0; i < len; i++) {
    switch ((char) block->seq.s[i]) {
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
            " A, T, C, G, or N: '%c', from %s", block->seq.s[i], block->seq.s);
    }
  }
}

void update_qual_matrices(kseq_t *block, int *qual_matrix, quality_type q_type) {
  /*
    Given `fastx_block` (must be type FASTQ), adjust the quality
    frequency in `qual_matrix` accordingly.
  */
  int i, len;
  int q_range = quality_contants[q_type][Q_MAX] - quality_contants[q_type][Q_MIN];
  int q_min = quality_contants[q_type][Q_MIN];
  int q_max = quality_contants[q_type][Q_MAX];
  int q_offset = quality_contants[q_type][Q_OFFSET];

  if (!block->qual.s)
    error("update_qual_matrices only works on FASTQ files");
  
  len = strlen(block->seq.s);
  if (len != strlen(block->qual.s))
    error("improperly formatted FASTQ file; sequence and quality lengths differ");

  for (i = 0; i < len; i++) {
    if ((char) block->qual.s[i] - q_offset < q_min || (char) block->qual.s[i] - q_offset > q_max)
      error("base quality out of range (%d < b < %d) encountered: %d", q_min,
            q_max, (char) block->qual.s[i]);
    
    qual_matrix[(q_range+1)*i + ((char) block->qual.s[i]) - q_offset - q_min]++;
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
  int is_missing, ret, size_out_list = 3, l;
  unsigned int num_unique_seqs = 0;

  unsigned int nblock = 0;
  kseq_t *block;
  SEXP base_counts, qual_counts, seq_hash, seq_hash_names, out_list, seq_lengths;
  int *ibc, *iqc, *isl, i, j, q_type, q_range;

  q_type = INTEGER(quality_type)[0];
  q_range = quality_contants[q_type][Q_MAX] - quality_contants[q_type][Q_MIN];
  
  if (LOGICAL(hash)[0]) {
    h = kh_init(str);
    kh_resize(str, h, 1572869);
    size_out_list = 4;
  }

  gzFile *fp = gzopen(CHAR(STRING_ELT(filename, 0)), "r");
  if (fp == NULL)
    error("failed to open file '%s'", CHAR(STRING_ELT(filename, 0)));
  block = kseq_init(fp);

  PROTECT(out_list = allocVector(VECSXP, size_out_list));
  PROTECT(base_counts = allocMatrix(INTSXP, NUM_BASES, INTEGER(max_length)[0]));
  PROTECT(qual_counts = allocMatrix(INTSXP, q_range + 1, INTEGER(max_length)[0]));
  PROTECT(seq_lengths = allocVector(INTSXP, INTEGER(max_length)[0]));
  
  ibc = INTEGER(base_counts);
  iqc = INTEGER(qual_counts);
  isl = INTEGER(seq_lengths);

  zero_int_matrix(ibc, NUM_BASES, INTEGER(max_length)[0]);
  zero_int_matrix(iqc, q_range + 1, INTEGER(max_length)[0]);

  while ((l = kseq_read(block)) >= 0) {
    R_CheckUserInterrupt();
    
    if (l >= LINE_BUFFER-1)
      error("read in sequence greater than LINE_BUFFER size");

    printf("name: %s\n", block->name.s); 
    if (block->comment.l) printf("comment: %s\n", block->comment.s);
    printf("seq: %s\n", block->seq.s);
    if (block->qual.l) printf("qual: %s\n", block->qual.s);
    
    update_base_matrices(block, ibc);
    update_qual_matrices(block, iqc, q_type);
    
    isl[nblock] = strlen(block->seq.s); //TODO use kseq's length

    if (LOGICAL(hash)[0]) {
      k = kh_get(str, h, block->seq.s);
      is_missing = (k == kh_end(h));
      if (is_missing) {
        k = kh_put(str, h, strdup(block->seq.s), &ret);
        kh_value(h, k) = 1;
        num_unique_seqs++;
      } else
        kh_value(h, k) = kh_value(h, k) + 1;
      
      if (LOGICAL(verbose)[0] && nblock % 100000 == 0)
        printf("on block %d, %d entries in hash table\n", nblock, num_unique_seqs);
    }
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
  block = kseq_init(fp);
  gzclose(fp);

  return out_list;
}
