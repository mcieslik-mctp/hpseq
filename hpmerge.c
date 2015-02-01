#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <zlib.h>
#include <kseq.h>

#define MAX_READ_LENGTH 10000

KSEQ_INIT(gzFile, gzread)


int main(int argc, char *argv[]) {
  // aux variables
  int c, i, j, k;
  int status = 0;
  char *fo = NULL;
  int drop = 0;
  int revc = 1;
  
  while ((c = getopt(argc, argv, "o:nd")) >= 0) {
    if (c == 'o') {
      fo = malloc((strlen(optarg) + 1));
      if (fo) strcpy(fo, optarg);
    }
    if (c == 'n') {
      revc = 0;
    }
    if (c == 'd') {
      drop = 1;
    }
  }

  if ((argv[optind] == NULL) || (argv[optind + 1] == NULL) || (fo == NULL)) {
  	fprintf(stderr, "\n");
  	fprintf(stderr,
  	"Usage: hpmerge -o <output.fq> <input_1.fq> <input_2.fq>\n\n");
  	fprintf(stderr, "Options:\n");
  	fprintf(stderr, "    -o STR    output file name\n");
  	fprintf(stderr, "    -o STR    output file name\n");
  	fprintf(stderr, "    -d        drop quality scores\n");
  	fprintf(stderr, "    -n        do not reverse complement\n");
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
    int last;
    FILE *fw;
    gzFile  fp[2];
    kseq_t *ks[2];
    char *fs; // forward sequence
    char rc[MAX_READ_LENGTH] = {0}; // reverse-complement
    
    fw = fopen(fo, "wb");
    fp[0] = gzopen(argv[optind], "rb");
    fp[1] = gzopen(argv[optind+1], "rb");
    ks[0] = kseq_init(fp[0]);
    ks[1] = kseq_init(fp[1]);
    while ((kseq_read(ks[0]) >= 0) && (kseq_read(ks[1]) >= 0)) {
      fputc((!drop) && ks[0]->qual.l? '@' : '>', fw);
      fputs(ks[0]->name.s, fw);
      fputc('\n', fw);
      fputs(ks[0]->seq.s, fw);
      // the second read is reverse-complemented
      if (revc) {
	last = ks[1]->seq.l-1;
	fs = ks[1]->seq.s;
	for (i=0; i<=last; ++i) {
	  switch(fs[i]) {
	  case 'A': rc[last - i]='T'; break;
	  case 'T': rc[last - i]='A'; break;
	  case 'G': rc[last - i]='C'; break;
	  case 'C': rc[last - i]='G'; break;
	  default:  rc[last - i]='N'; break;
	  }
	}
	fputs(rc, fw);
      } else {
	fputs(ks[1]->seq.s, fw);
      }
      fputc('\n', fw);
      if ((!drop) && ks[0]->qual.l) {
	fputs("+", fw);
	fputc('\n', fw);
	fputs(ks[0]->qual.s, fw);
	fputs(ks[1]->qual.s, fw);
	fputc('\n', fw);
      }
    }
    // cleanup
    kseq_destroy(ks[0]);
    kseq_destroy(ks[1]);
    gzclose(fp[0]);
    gzclose(fp[1]);
    fclose(fw);
  }

  /* // cleanup */
  free(fo);
  return status;
}
