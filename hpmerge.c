#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <zlib.h>
#include <kseq.h>

#define MAX_READ_LENGTH 10000

KSEQ_INIT(gzFile, gzread)

int HammingDistance(const char* a, int na, int mm, const char* b) {
  int num_mismatches = 0;
  while ((na > 0) && (num_mismatches <= mm)) {
    if (*a != *b)
      ++num_mismatches;
    --na;
    ++a;
    ++b;
  }
  return num_mismatches;
}

int main(int argc, char *argv[]) {
  // aux variables
  int status = 0;
  int c, i, j, k;
  char buf[10000];
  // malloc'ed
  char *fo[3] = {NULL, NULL, NULL};

  // parameters
  int p_drop = 0; // drop quality scores
  int p_min_off = 10; // minimum offset == minimum overlap - 1
  double p_max_err = 0.12; // maximum error

  while ((c = getopt(argc, argv, "o:m:e:d")) >= 0) {
    if (c == 'o') {
      snprintf(buf, sizeof(buf), "%s%s", optarg, "_merged.fq");
      fo[0] = malloc((strlen(buf) + 1));
      if (fo[0]) strcpy(fo[0], buf);
      snprintf(buf, sizeof(buf), "%s%s", optarg, "_unmerged_1.fq");
      fo[1] = malloc((strlen(buf) + 1));
      if (fo[1]) strcpy(fo[1], buf);
      snprintf(buf, sizeof(buf), "%s%s", optarg, "_unmerged_2.fq");
      fo[2] = malloc((strlen(buf) + 1));
      if (fo[2]) strcpy(fo[2], buf);
    }
    else if (c == 'm') p_min_off = atoi(optarg) - 1;
    else if (c == 'e') p_max_err = atof(optarg);
    else if (c == 'd') p_drop = 1;
  }
  if ((argv[optind] == NULL) || (argv[optind + 1] == NULL) || (fo[0] == NULL)) {
    fprintf(stderr, "\n");
    fprintf(stderr,
    "Usage: hpmerge [-d] [-e] [-m] -o pfx <input_1.fq> <input_2.fq>\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o STR    output file prefix\n");
    fprintf(stderr, "    -m INT    minimum overlap\n");
    fprintf(stderr, "    -e FLOAT  maximum error\n");
    fprintf(stderr, "    -d        drop quality scores\n");
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
    int mm;
    int hd;
    int max_over; // maximum overlap
    int last0, last1; // index of last base
    char *s0, *s1;  // sequence
    char *c1, *rc1; // complement and reverse-complement
    c1  = calloc(1000 * sizeof(char), sizeof(char));
    rc1 = calloc(1000 * sizeof(char), sizeof(char));
    FILE *fw[3];
    gzFile  fp[2];
    kseq_t *ks[2];
    kseq_t *r;
    FILE *f;
    int merged;
    // file io
    fw[0] = fopen(fo[0], "wb");
    fw[1] = fopen(fo[1], "wb");
    fw[2] = fopen(fo[2], "wb");
    fp[0] = gzopen(argv[optind], "rb");
    fp[1] = gzopen(argv[optind+1], "rb");
    ks[0] = kseq_init(fp[0]);
    ks[1] = kseq_init(fp[1]);
    while ((kseq_read(ks[0]) >= 0) && (kseq_read(ks[1]) >= 0)) {
      merged = 0;
      //
      last1 = ks[1]->seq.l-1;
      s1 = ks[1]->seq.s;
      // complement and reverse complement
      for (i=0; i<=last1; ++i) {
        switch(s1[i]) {
        case 'A': c1[i] = rc1[last1 - i] = 'T'; break;
        case 'T': c1[i] = rc1[last1 - i] = 'A'; break;
        case 'G': c1[i] = rc1[last1 - i] = 'C'; break;
        case 'C': c1[i] = rc1[last1 - i] = 'G'; break;
        default:  c1[i] = rc1[last1 - i] = 'N'; break;
        }
      }
      s1 = rc1;
      c1[i+1] = '\0';
      last0 = ks[0]->seq.l-1;
      s0 = ks[0]->seq.s;
      //
      int offset = (last0 < last1 ? last0 : last1);
      for (i=offset; i>=p_min_off; --i) {
        mm = (int)((i+1) * p_max_err);
        hd = HammingDistance(&s0[last0 - i], i+1, mm, s1);
        f = fw[0];
        if (hd <= mm) {
          // merge
          fputc((!p_drop) && ks[0]->qual.l? '@' : '>', fw[0]);
          fputs(ks[0]->name.s, fw[0]);
          fputc('\n', fw[0]);
          s0[last0 - i] = '\0';
          fputs(s0, fw[0]);
          fputs(&(c1[i + 1]), fw[0]);
          fputc('\n', fw[0]);
          if ((!p_drop) && ks[0]->qual.l) {
            fputs("+", fw[0]);
            fputc('\n', fw[0]);
            ks[0]->qual.s[last0 - i] = '\0';
            fputs(ks[0]->qual.s, fw[0]);
            fputs(&ks[1]->qual.s[i+1], fw[0]);
            fputc('\n', fw[0]);
          }
          merged = 1;
          break;
        }
      }
      if (!merged) {
  	for (i = 1; i < 3; ++i) {
  	  f = fw[i];
  	  r = ks[i-1];
  	  fputc((!p_drop) && r->qual.l? '@' : '>', f);
  	  fputs(r->name.s, f);
  	  fputc('\n', f);
  	  fputs(r->seq.s, f);
  	  fputc('\n', f);
  	  if ((!p_drop) && r->qual.l) {
  	    fputs("+", f);
  	    fputc('\n', f);
  	    fputs(r->qual.s, f);
  	    fputc('\n', f);
  	  }
  	}
      }
    }
    // cleanup
    free(c1);
    free(rc1);
    kseq_destroy(ks[0]);
    kseq_destroy(ks[1]);
    gzclose(fp[0]);
    gzclose(fp[1]);
    fclose(fw[0]);
    fclose(fw[1]);
    fclose(fw[2]);
    
  }
  // cleanup
  free(fo[0]);
  free(fo[1]);
  free(fo[2]);
  return status;
}
