#include <config.h>
#include <stdio.h>
/* #include <string.h> */
/* #include <stdlib.h> */
/* #include <unistd.h> */
#include <getopt.h>
#include <stdint.h>
#include <zlib.h>
#include <kwset.h>
#include <kseq.h>

KSEQ_INIT(gzFile, gzread)


int main(int argc, char *argv[]) {
  // aux variables
  int c, i, j, k;
  int status = 0;
  char* fo[2] = {NULL, NULL};
    
  kwset_t seqs = kwsalloc(NULL);
  int nseqs = 0;
  while ((c = getopt(argc, argv, "1:2:s:")) >= 0) {
    if (c == 's') {
      nseqs++;
      kwsincr(seqs, optarg, strlen(optarg));
    }
    if (c == '1') {
      fo[0] = malloc((strlen(optarg) + 1));
      if (fo[0]) strcpy(fo[0], optarg);
    }
    if (c == '2') {
      fo[1] = malloc((strlen(optarg) + 1));
      if (fo[1]) strcpy(fo[1], optarg);
    }
  }
  kwsprep(seqs);

  if ((argv[optind] == NULL) || (argv[optind + 1] == NULL) ||
      (fo[0] == NULL) || (fo[1] == NULL) || (nseqs == 0)
      ) {
  	fprintf(stderr, "\n");
  	fprintf(stderr,
  	"Usage: hpscan_cw [options] -s <query> [-s <query>] -1 <o_1.fq> -2 <o_2.fq> <i_1.fq> <i_2.fq>\n\n");
  	fprintf(stderr, "Options:\n");
  	fprintf(stderr, "    -s STR    query sequence(s)\n");
  	fprintf(stderr, "    -1 STR    output file name for read 1\n");	
  	fprintf(stderr, "    -2 STR    output file name for read 2\n");
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
    gzFile  fp[2];
    kseq_t *ks[2];
    fp[0] = gzopen(argv[optind], "rb");
    fp[1] = gzopen(argv[optind+1], "rb");
    ks[0] = kseq_init(fp[0]);
    ks[1] = kseq_init(fp[1]);
    struct kwsmatch kwsm;
    int match;
    long res;
    kseq_t *r;
    FILE *f;
    while ((kseq_read(ks[0]) >= 0) && (kseq_read(ks[1]) >= 0)) {
      match = 0;
      for (i = 0; i < 2 && !match; ++i) {
	res = (long)kwsexec(seqs, ks[i]->seq.s, ks[i]->seq.l, &kwsm);
	if (res >= 0) {
	  match = 1;
	  break;
	}
      }
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
    kseq_destroy(ks[0]);
    kseq_destroy(ks[1]);
    gzclose(fp[0]);
    gzclose(fp[1]);
    fclose(fw[0]);
    fclose(fw[1]);
  }

  // cleanup
  kwsfree(seqs);
  if (fo[0] != NULL) {
    free(fo[0]);
  }
  if (fo[1] != NULL) {
    free(fo[1]);
  }
  
  return status;
}
