#include "demultiplex.h"
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

void opt_init(opt_t *opt)
{
  memset(opt, 0, sizeof(opt_t));

  opt->LEN = 150; opt->MASK_LEN = opt->LEN / 2;

  opt->match = 2; opt->mismatch = 2;
  opt->open = 3; opt->ext = 1;

  opt->n_threads = 3;
  opt->chunk_size = 400000000;

  opt->compress = 0;
  opt->output_unclassified = 1;

  opt->bc = (bc_t *)malloc(sizeof(bc_t));
  memset(opt->bc, 0, sizeof(bc_t));
}

void opt_set_mat(int match, int mismatch, int8_t mat[25])
{
  int i, j, k;
  for (i = k = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      mat[k++] = i == j? match : -mismatch;
    }
    mat[k++] = 0;
  }
  for (j = 0; j < 5; ++j) mat[k++] = 0;
}

void opt_allocate_bc(opt_t *opt, int8_t *len)
{
  bc_t *bc = opt->bc;
  bc->capacity = 16;
  bc->name = (char **)malloc(sizeof(char *) * 16);
  bc->len = (int8_t *)malloc(sizeof(int8_t) * 16);
  bc->nt = (int8_t **)malloc(sizeof(int8_t *) * 16);
  bc->nt_rc = (int8_t **)malloc(sizeof(int8_t *) * 16);
  gzFile fp = gzopen(opt->fb, "r");
  kseq_t *seq = kseq_init(fp);
  while(kseq_read(seq) >= 0) {
    if (bc->idx + 1 >= bc->capacity) {
      bc->capacity <<= 1;
      bc->name = (char **)realloc(bc->name, sizeof(char *) * bc->capacity);
      bc->len = (int8_t *)realloc(bc->len, sizeof(int8_t) * bc->capacity);
      bc->nt = (int8_t **)realloc(bc->nt, sizeof(int8_t *) * bc->capacity);
      bc->nt_rc = (int8_t **)realloc(bc->nt_rc, sizeof(int8_t *) * bc->capacity);
    }
    if (seq->seq.l < *len) *len = seq->seq.l;
    bc->name[bc->idx] = strdup(seq->name.s);
    bc->len[bc->idx] = seq->seq.l;
    bc->nt[bc->idx] = seq_to_nt(seq->seq.s, seq->seq.l, 0);
    bc->nt_rc[bc->idx] = seq_to_nt(seq->seq.s, seq->seq.l, 1);
    bc->idx += 1;
  }
  kseq_destroy(seq);
  gzclose(fp);

  bc->file_num = bc->idx;
}

void opt_hash_bc(opt_t *opt, char **name)
{
  int i = 0;

  khint_t k; int absent;
  opt->bc->h = kh_init(dual);

  FILE *fp = fopen(opt->dual, "r");
  char *line = NULL; size_t len;
  char fn[128], bc1[128], bc2[128], tmp[256];
  while(getline(&line, &len, fp) != -1) {
    sscanf(line, "%s\t%s\t%s", fn, bc1, bc2);
    sprintf(tmp, "%s%s", bc1, bc2);
    tmp[strlen(bc1) + strlen(bc2)] = '\0';
    k = kh_put(dual, opt->bc->h, tmp, &absent);
    if (absent) {
      name[i] = strdup(fn);
      kh_key(opt->bc->h, k) = strdup(tmp);
      kh_val(opt->bc->h, k) = i++;
    }
  }
  free(line);
  fclose(fp);

  opt->bc->file_num = i;
}

void opt_set_output(opt_t *opt, char **name)
{
  if (access(opt->path, F_OK)) {
    int status = mkdir(opt->path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status) print_log("Failed to create output dir %s.", opt->path), exit(1);
  }

  int i; char tmp[1024];

  bc_t *bc = opt->bc;
  bc->compressed = opt->compress;
  bc->total_files = bc->file_num + (opt->output_unclassified ? 1 : 0);
  if (bc->compressed) {
    bc->gptr = (gzFile *)malloc(sizeof(gzFile) * bc->total_files);
    bc->ptr = NULL;
  } else {
    bc->ptr = (FILE **)malloc(sizeof(FILE *) * bc->total_files);
    bc->gptr = NULL;
  }
  bc->buffer = (char **)malloc(sizeof(char *) * bc->total_files);
  bc->offset = (int32_t *)calloc(bc->total_files, sizeof(int32_t));

  const char *ext = opt->mode ? ".fasta" : ".fastq";
  const char *suffix = opt->compress ? ".gz" : "";

  for(i = 0; i < bc->file_num; i++) {
    snprintf(tmp, sizeof(tmp), "%s/%s%s%s", opt->path, name[i], ext, suffix);
    if (bc->compressed) {
      bc->gptr[i] = gzopen(tmp, "wb");
      if (!bc->gptr[i]) {
        print_log("Failed to open output file %s: %s", tmp, strerror(errno));
        exit(1);
      }
    } else {
      bc->ptr[i] = fopen(tmp, "w");
      if (!bc->ptr[i]) {
        print_log("Failed to open output file %s: %s", tmp, strerror(errno));
        exit(1);
      }
    }
    bc->buffer[i] = (char *)malloc(sizeof(char) * BUFSIZE);
  }

  if (opt->output_unclassified) {
    snprintf(tmp, sizeof(tmp), "%s/unclassified%s%s", opt->path, ext, suffix);
    int idx = bc->file_num;
    if (bc->compressed) {
      bc->gptr[idx] = gzopen(tmp, "wb");
      if (!bc->gptr[idx]) {
        print_log("Failed to open output file %s: %s", tmp, strerror(errno));
        exit(1);
      }
    } else {
      bc->ptr[idx] = fopen(tmp, "w");
      if (!bc->ptr[idx]) {
        print_log("Failed to open output file %s: %s", tmp, strerror(errno));
        exit(1);
      }
    }
    bc->buffer[idx] = (char *)malloc(sizeof(char) * BUFSIZE);
  }
}

void opt_set_bc(opt_t *opt)
{
  int i;
  int8_t len = INT8_MAX;

  opt_allocate_bc(opt, &len);

  if (opt->dual) {
    char **name = (char **)malloc(sizeof(char *) * 1024);

    opt_hash_bc(opt, name);
    opt_set_output(opt, name);

    for(i = 0; i < opt->bc->file_num; i++)  free(name[i]);
    free(name);
  } else opt_set_output(opt, opt->bc->name);

  if (!opt->flag) {
    int gap = (int)(len * ERR + 1);
    int match = len - gap;
    opt->score = match * opt->match - gap * opt->open;
  }
  print_log("Minimal alignment score is %d", opt->score);
}

void opt_set(opt_t *opt, char *fn)
{
  opt_set_mat(opt->match, opt->mismatch, opt->mat);

  gzFile fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
  opt->ks = kseq_init(fp);

  opt_set_bc(opt);
}

void opt_free(opt_t *opt)
{
  int i;
  bc_t *bc = opt->bc;

  for(i = 0; i < bc->total_files; i++) {
    bc_flush(opt, i);
    if (bc->compressed) gzclose(bc->gptr[i]);
    else fclose(bc->ptr[i]);
    free(bc->buffer[i]);
  }
  for(i = 0; i < bc->idx; i++) {
    free(bc->name[i]);
    free(bc->nt[i]); free(bc->nt_rc[i]);
  }
  if (opt->dual) {
    khint_t k;
    for(k = kh_begin(bc->h); k != kh_end(bc->h); ++k)
      if (kh_exist(bc->h, k)) free((char*)kh_key(bc->h, k));
    kh_destroy(dual, bc->h);
    free(opt->dual);
  }
  if (bc->ptr) free(bc->ptr);
  if (bc->gptr) free(bc->gptr);
  free(bc->buffer); free(bc->offset);
  free(bc->name); free(bc->len);
  free(bc->nt); free(bc->nt_rc);
  free(bc);

  gzFile fp = opt->ks->f->f;
  kseq_destroy(opt->ks);
  gzclose(fp);
  free(opt->fb); free(opt->path);
  if (opt->log) fclose(opt->log);
}
