#include <stdio.h>
/* #include <string.h> */
/* #include <stdlib.h> */
/* #include <unistd.h> */
#include <getopt.h>
/* #include <stdint.h> */
#include <zlib.h>
#include <kseq.h>
#include <ksw.h>
KSEQ_INIT(gzFile, gzread)

#define MAX_SEQS 100

// letters to numbers
unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct ta_seq_t{
	size_t len;
	uint8_t *seq;
        kswq_t *qp;
	uint64_t cnt;
} ta_seq_t;


int main(int argc, char *argv[]) {
  // aux variables
  int c, i, j, k;
  int status = 0;

  // alignment score defaults
  int ma = 1, mm = 2, go = 0, ge = 2;
  int min_score = 13;

  //
  char* fo[2] = {NULL, NULL};
    
  ta_seq_t *seqs;
  seqs = malloc(MAX_SEQS * sizeof(ta_seq_t)); 

  int nseqs = 0;
  while ((c = getopt(argc, argv, "1:2:s:m:n:o:e:")) >= 0) {
    if (c == 's') {
      ta_seq_t *ps;
      ps = &seqs[nseqs++];
      char* pc = malloc((strlen(optarg) + 1) * sizeof(char));
      if (pc) strcpy(pc, optarg);
      ps->seq = (uint8_t*) pc;
      ps->qp = 0;
      ps->len = (int)strlen(pc);
      for (int i = 0; i < ps->len; ++i)
	ps->seq[i] = seq_nt4_table[(uint8_t)ps->seq[i]];
    }
    if (c == '1') {
      fo[0] = malloc((strlen(optarg) + 1) * sizeof(char));
      if (fo[0]) strcpy(fo[0], optarg);
    }
    if (c == '2') {
      fo[1] = malloc((strlen(optarg) + 1) * sizeof(char));
      if (fo[1]) strcpy(fo[1], optarg);
    }
    switch (c) {
    case 'a': ma = atoi(optarg); break;
    case 'm': mm = atoi(optarg); break;
    case 'o': go = atoi(optarg); break;
    case 'e': ge = atoi(optarg); break;
    case 'l': ge = atoi(optarg); break;
    }
  }

  /* int got_pipe; */
  /* got_pipe = !isatty(STDIN_FILENO); */
  if ((argv[optind] == NULL) || (argv[optind + 1] == NULL) ||
      (fo[0] == NULL) || (fo[1] == NULL) || (nseqs == 0)
      ) {
  	fprintf(stderr, "\n");
  	fprintf(stderr, "Usage:  swfq [options] -s <query> -1 <o_1.fq> -2 <o_2.fq> <i_1.fq> <i_2.fq>\n\n");
  	fprintf(stderr, "Options:\n");
  	fprintf(stderr, "    -s STR    query sequence(s)\n");
  	fprintf(stderr, "    -l INT    min length (score): %d\n", min_score);
  	fprintf(stderr, "    -a INT    match score: %d\n", ma);
  	fprintf(stderr, "    -m INT    mis-match penalty: %d\n", mm);
  	fprintf(stderr, "    -o INT    gap-open penalty: %d\n", go);
  	fprintf(stderr, "    -e INT    gap-extend penalty: %d\n", ge);
  	fprintf(stderr, "\n");
	status = 1;
  } else {
    if(access(argv[optind], F_OK) == -1) {
      printf("missing file %s\n", argv[optind]);
      status = 1;
    }
    if(access(argv[optind+1], F_OK) == -1) {
      printf("missing file %s\n", argv[optind+1]);
      status = 1;
    }
  }
  

  if (!status) {
    FILE *fw[2];
    fw[0] = fopen(fo[0], "wb");
    fw[1] = fopen(fo[1], "wb");
    
    gzFile fp[2];
    kseq_t *ks[2];
    fp[0] = gzopen(argv[optind], "rb");
    fp[1] = gzopen(argv[optind+1], "rb");
    ks[0] = kseq_init(fp[0]);
    ks[1] = kseq_init(fp[1]);
    
    // DP score-matrix
    int8_t mat[25];
    for (i = k = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j)
      mat[k++] = (i == j) ? ma : -mm;
      mat[k++] = 0; // ambiguous base
    }
    for (j = 0; j < 5; ++j) mat[k++] = 0;
    //
    kstring_t str = {0,0,0};
    while ((kseq_read(ks[0]) >= 0) && (kseq_read(ks[1]) >= 0)) {
      int match = 0;
      for (i = 0; i < 2 && !match; ++i) {
	if (str.m < ks[i]->seq.m) {
	  str.m = ks[i]->seq.m;
	  str.s = realloc(str.s, str.m);
	}
	str.l = ks[i]->seq.l;
	for (j = 0; j < ks[i]->seq.l; ++j) {
	  str.s[j] = seq_nt4_table[(uint8_t)ks[i]->seq.s[j]];
	}
	for (j = 0; j < nseqs && !match; ++j) {
	  kswr_t swa;
	  ta_seq_t *ps = &seqs[j];
	  swa = ksw_align(ps->len, ps->seq, str.l, (uint8_t*)str.s, 5, mat, go, ge,
			  KSW_XBYTE|KSW_XSTART|min_score, &ps->qp);
	  if (swa.score > min_score) {
	    match = 1;
	    break;
	  }
	}
      }
      // 
      kseq_t *r;
      FILE *f;
      if (match) {
	for (i = 0; i < 2; ++i) {
	  f = fw[i];
	  r = ks[i];
	  fputc(r->qual.l? '@' : '>', f);
	  fputs(r->name.s, f);
	  fputc('\n', f);
	  fputs(r->seq.s, f);
	  fputc('\n', f);
	  if (r->qual.l) {
	    fputs("+", f);
	    fputc('\n', f);
	    fputs(r->qual.s, f);
	    fputc('\n', f);
	  }
	}
      }
    }
    
    // cleanup
    free(str.s);
    kseq_destroy(ks[0]);
    kseq_destroy(ks[1]);
    gzclose(fp[0]);
    gzclose(fp[1]);
    fclose(fw[0]);
    fclose(fw[1]);
  }

  // cleanup
  for (int j = 0; j < nseqs; ++j) {
    ta_seq_t *p = &seqs[j];
    free(p->seq);
    free(p->qp);
  }
  free(seqs);
  if (fo[0] != NULL) {
    free(fo[0]);
  }
  if (fo[1] != NULL) {
    free(fo[1]);
  }
  return status;
}
