#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <zlib.h>
#include <kseq.h>

KSEQ_INIT(gzFile, gzread)


int main(int argc, char *argv[]) {
  // aux variables
  int c, i, j, k;
  int status = 0;
  char *fo = NULL;
  int drop = 0;
  
  while ((c = getopt(argc, argv, "o:d")) >= 0) {
    if (c == 'o') {
      fo = malloc((strlen(optarg) + 1));
      if (fo) strcpy(fo, optarg);
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
    FILE *fw;
    gzFile  fp[2];
    kseq_t *ks[2];
    
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
      fputs(ks[1]->seq.s, fw);
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
