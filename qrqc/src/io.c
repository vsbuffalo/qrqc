#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include "samtools/khash.h"
#include "samtools/kseq.h"
#include "io.h"

#ifdef _WIN32
int gzreadclone(FILE *file, void *buf, unsigned int len) {
  unsigned int i;
  char c;
  for (i = 0; i < len; i++) {
    c = getc(file);
    if (c == EOF) return i;
    ((char *)buf)[i] = c;
  }
  return i;
}
#define FILE_TYPE FILE
#define FILE_OPEN(x, m) (fopen(x, m))
#define FILE_CLOSE(x) (fclose(x))
KSEQ_INIT(FILE_TYPE*, gzreadclone)
#else
#include <zlib.h>
#define FILE_TYPE gzFile
#define FILE_OPEN(x, m) (gzopen(x, m))
#define FILE_CLOSE(x) (gzclose(x))
KSEQ_INIT(gzFile, gzread)
#endif

KHASH_MAP_INIT_STR(str, int)

#define INIT_MAX_SEQ 500
#define NUM_BASES 16 /* includes all IUPAC codes. */

#define IS_FASTQ(quality_type) INTEGER(quality_type)[0] >= 0

typedef enum {
  PHRED,
  SANGER,
  SOLEXA,
  ILLUMINA
} quality_type;

#define Q_OFFSET 0
#define Q_MIN 1
#define Q_MAX 2

static const int quality_contants[4][3] = {
  /* offset, min, max */
  {0, 4, 60}, /* PHRED */
  {33, 0, 93}, /* SANGER */
  {64, -5, 62}, /* SOLEXA */
  {64, 0, 62} /* ILLUMINA */
};


static void update_base_matrices(kseq_t *block, int *base_matrix) {
  /*
    Given a *kseq_t block containing a sequence, update the
    *base_matrix to contain this block's bases (by position).
  */
  int i;
  
  for (i = 0; i < block->seq.l; i++) {
    switch ((char) block->seq.s[i]) {
    case 'A':
      base_matrix[NUM_BASES*i]++;
      break;
    case 'T':
      base_matrix[NUM_BASES*i + 1]++;
      break;
    case 'C':
      base_matrix[NUM_BASES*i + 2]++;
      break;
    case 'G':
      base_matrix[NUM_BASES*i + 3]++;
      break;
    case 'N':
      base_matrix[NUM_BASES*i + 4]++;
      break;
    case 'R':
      base_matrix[NUM_BASES*i + 5]++;
      break;
    case 'Y':
      base_matrix[NUM_BASES*i + 6]++;
      break;
    case 'S':
      base_matrix[NUM_BASES*i + 7]++;
      break;
    case 'W':
      base_matrix[NUM_BASES*i + 8]++;
      break;
    case 'K':
      base_matrix[NUM_BASES*i + 9]++;
      break;
    case 'M':
      base_matrix[NUM_BASES*i + 10]++;
      break;
    case 'B':
      base_matrix[NUM_BASES*i + 11]++;
      break;
    case 'D':
      base_matrix[NUM_BASES*i + 12]++;
      break;
    case 'H':
      base_matrix[NUM_BASES*i + 13]++;
      break;
    case 'V':
      base_matrix[NUM_BASES*i + 14]++;
      break;
    case '.':
    case '-':
      base_matrix[NUM_BASES*i + 15]++;
      break;
    default:
      error("Sequence character encountered that is not"
            "a IUPAC code: '%c', from %s", block->seq.s[i], block->seq.s);
    }
  }
}

static void update_qual_matrices(kseq_t *block, int *qual_matrix, quality_type q_type) {
  /*
    Given a *kseq_t block containing qualities, update the
    *qual_matrix to contain this block's qualities (by position).
  */
  int i;
  int q_range = quality_contants[q_type][Q_MAX] - quality_contants[q_type][Q_MIN];
  int q_min = quality_contants[q_type][Q_MIN];
  int q_max = quality_contants[q_type][Q_MAX];
  int q_offset = quality_contants[q_type][Q_OFFSET];

  if (!block->qual.l && block->seq.l > 0)
    error("update_qual_matrices only works on FASTQ files");
  
  for (i = 0; i < block->qual.l; i++) {
    R_CheckUserInterrupt();
    if ((char) block->qual.s[i] - q_offset <= q_min || (char) block->qual.s[i] - q_offset >= q_max)
      error("base quality out of range (%d <= b <= %d) encountered: %d", q_min,
            q_max, (char) block->qual.s[i]);
    
    qual_matrix[(q_range+1)*i + ((char) block->qual.s[i]) - q_offset - q_min]++;
  }
}

static void zero_int_matrix(int *matrix, int nx, int ny) {
  /*
    Zero out a matrix of integers.
  */
  int i, j;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++)
      matrix[i + nx*j] = 0;
  }
}

static void zero_int_vector(int *vector, int n) {
  /*
    Zero out a vector of integers.
  */
  int i;
  for (i = 0; i < n; i++) {
    vector[i] = 0;
  }
}

static void add_seq_to_khash(khash_t(str) *h, kseq_t *block, unsigned int *num_unique_seqs) {
  /*
    Given a hash of strings (khash_t(str) *h), a block containing a
    sequence, update the hash. If the sequence exist, increase its
    value by 1. If not, add the sequence to the hash and increase
    num_unique_seqs.
  */
  khiter_t k;
  int is_missing, ret;
  k = kh_get(str, h, block->seq.s);
  is_missing = (k == kh_end(h));
  if (is_missing) {
    k = kh_put(str, h, strdup(block->seq.s), &ret);
    kh_value(h, k) = 1;
    (*num_unique_seqs)++;
  } else
    kh_value(h, k) = kh_value(h, k) + 1;
}

static void seq_khash_to_VECSXP(khash_t(str) *h, SEXP seq_hash, SEXP seq_hash_names) {
  /*
    Given a hash with string keys, output values to a given
    (pre-allocated!) SEXP seq_hash. Then, output keys to
    seq_hash_names.
   */
  khiter_t k;
  int i;
  
  i = 0;
  for (k = kh_begin(h); k != kh_end(h); ++k) {
    R_CheckUserInterrupt();
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
}

extern SEXP summarize_file(SEXP filename, SEXP max_length, SEXP quality_type, SEXP hash, SEXP verbose) {
  /*
    Given a FASTA or FASTQ file, read in sequences and gather
    statistics on bases, qualities, sequence lengths, and unique
    sequences (if hash is true). 

    All matrices are pre-allocated to max_length, and then trimmed
    accordingly in R. This may be (and should be) changed in future
    versions.

    Note that quality_type is -1 if the file is FASTA.
   */
  if (!isString(filename))
    error("filename should be a string");

  khash_t(str) *h=NULL;
  kseq_t *block;
  int size_out_list = 4, l, protects=0;
  unsigned int num_unique_seqs = 0, nblock = 0;
  /* Note: NULL and 0 initializations to stop warnings on Windows systems */
  SEXP base_counts, out_list, seq_lengths, qual_counts=NULL, seq_hash=NULL, seq_hash_names=NULL;
  int *ibc, *isl, *iqc=NULL, q_type=0, q_range=0;

  if (IS_FASTQ(quality_type)) {
    q_type = INTEGER(quality_type)[0];
    q_range = quality_contants[q_type][Q_MAX] - quality_contants[q_type][Q_MIN];
    size_out_list++;
  }
  
  if (LOGICAL(hash)[0]) {
    h = kh_init(str);
    kh_resize(str, h, 1572869);
    size_out_list++;
  }

  FILE_TYPE *fp = FILE_OPEN(CHAR(STRING_ELT(filename, 0)), "r");
  if (fp == NULL)
    error("failed to open file '%s'", CHAR(STRING_ELT(filename, 0)));
  block = kseq_init(fp);

  PROTECT(out_list = allocVector(VECSXP, size_out_list));
  PROTECT(base_counts = allocMatrix(INTSXP, NUM_BASES, INTEGER(max_length)[0]));
  PROTECT(seq_lengths = allocVector(INTSXP, INTEGER(max_length)[0]));
  protects = 3;

  ibc = INTEGER(base_counts);
  isl = INTEGER(seq_lengths);
  zero_int_matrix(ibc, NUM_BASES, INTEGER(max_length)[0]);
  zero_int_vector(isl, INTEGER(max_length)[0]);

  if (IS_FASTQ(quality_type)) {
    PROTECT(qual_counts = allocMatrix(INTSXP, q_range + 1, INTEGER(max_length)[0]));
    protects++;
    iqc = INTEGER(qual_counts);
    zero_int_matrix(iqc, q_range + 1, INTEGER(max_length)[0]);  
  }

  while ((l = kseq_read(block)) >= 0) {
    R_CheckUserInterrupt();
    if (IS_FASTQ(quality_type) && l == -2)
      error("improperly formatted FASTQ file; truncated quality string");
    if (l >= INTEGER(max_length)[0]-1)
      error("read in sequence greater than max.length");

    update_base_matrices(block, ibc);
    if (IS_FASTQ(quality_type))
      update_qual_matrices(block, iqc, q_type);
    
    isl[block->seq.l]++;

    if (LOGICAL(hash)[0]) {
      add_seq_to_khash(h, block, &num_unique_seqs);
      if (LOGICAL(verbose)[0] && nblock % 100000 == 0)
        printf("on block %d, %d entries in hash table\n", nblock, num_unique_seqs);
    }
    nblock++;
  }

  if (LOGICAL(hash)[0]) {
    PROTECT(seq_hash = allocVector(VECSXP, num_unique_seqs));
    PROTECT(seq_hash_names = allocVector(VECSXP, num_unique_seqs));
    protects += 2;
    if (LOGICAL(verbose)[0])
      printf("processing complete... now loading C hash structure to R...\n");
    seq_khash_to_VECSXP(h, seq_hash, seq_hash_names);
    kh_destroy(str, h);
  }

  SET_VECTOR_ELT(out_list, 0, base_counts);
  SET_VECTOR_ELT(out_list, 1, seq_lengths);
  if (IS_FASTQ(quality_type)) 
    SET_VECTOR_ELT(out_list, 2, qual_counts);
  
  if (LOGICAL(hash)[0]) {
    setAttrib(seq_hash, R_NamesSymbol, seq_hash_names);
    SET_VECTOR_ELT(out_list, 3, seq_hash);
  }

  /* One more protected SEXP from qual_counts */

  UNPROTECT(protects);

  block = kseq_init(fp);
  FILE_CLOSE(fp);

  return out_list;
}
