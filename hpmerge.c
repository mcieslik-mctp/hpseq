#include <stdio.h>
#include <getopt.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <kseq.h>

#define MAX_READ_LENGTH 1000
#define MAX_PATH_LENGTH 10000
#define PASS 0
#define CHOP 1
#define JOIN 2


KSEQ_INIT(gzFile, gzread)

inline int HammingDistance(const char* a, int na, int mm, const char* b) {
  int num_mismatches = 0;
  do {
    if (*(a++) != *(b++))
      ++num_mismatches;
    --na;
  } while ((num_mismatches <= mm) && (na > 0));
  return num_mismatches;
}

int main(int argc, char *argv[]) {
  // aux variables
  int status = 0;
  int c = 0, i = 0;
  int bt = 0;
  int verbose = 0;
  char buf[MAX_PATH_LENGTH];
  // malloc'ed
  char *fo[3] = {NULL, NULL, NULL};
  uint64_t opcounts[3] = {0};

  // parameters
  int p_chop = 1; // chop overlapping reads
  int p_join = 1; // join overlapping reads
  int p_min_off_chop =  5; // minimum offset == minimum overlap - 1
  double p_max_err_chop = 0.12; // maximum error
  int p_min_off_join = 14; // minimum offset == minimum overlap - 1
  double p_max_err_join = 0.12; // maximum error

  while ((c = getopt(argc, argv, "o:m:e:n:f:v:cj")) >= 0) {
    if (c == 'o') {
      snprintf(buf, sizeof(buf), "%s%s", optarg, "_0.fq");
      fo[0] = malloc((strlen(buf) + 1));
      if (fo[0]) strcpy(fo[0], buf);
      snprintf(buf, sizeof(buf), "%s%s", optarg, "_1.fq");
      fo[1] = malloc((strlen(buf) + 1));
      if (fo[1]) strcpy(fo[1], buf);
      snprintf(buf, sizeof(buf), "%s%s", optarg, "_2.fq");
      fo[2] = malloc((strlen(buf) + 1));
      if (fo[2]) strcpy(fo[2], buf);
    }
    else if (c == 'm') p_min_off_join = atoi(optarg) - 1;
    else if (c == 'e') p_max_err_join = atof(optarg);
    else if (c == 'n') p_min_off_chop = atoi(optarg) - 1;
    else if (c == 'f') p_max_err_chop = atof(optarg);
    else if (c == 'v') verbose = atoi(optarg);
    else if (c == 'j') p_join = 0;
    else if (c == 'c') p_chop = 0;
  }
  if ((argv[optind] == NULL) || (argv[optind + 1] == NULL) || (fo[0] == NULL) ||
      ((p_join && p_chop) && (p_max_err_chop < p_max_err_join)) || // chopping has to  
      ((p_join && p_chop) && (p_min_off_chop > p_min_off_chop))    // be more 'loose'
      ) {
    fprintf(stderr, "\n");
    fprintf(stderr,
    "Usage: hpmerge [-d] [-e] [-m] -o pfx <input_1.fq> <input_2.fq>\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o STR    output file prefix (required)\n");
    fprintf(stderr, "    -m INT    minimum join overlap [default: %d]\n", p_min_off_join + 1);
    fprintf(stderr, "    -e FLOAT  maximum join error [default: %.3f]\n", p_max_err_join);
    fprintf(stderr, "    -n INT    minimum chop overlap [default: %d]\n", p_min_off_chop + 1);
    fprintf(stderr, "    -f FLOAT  maximum chop error [default: %.3f]\n", p_max_err_chop);
    fprintf(stderr, "    -j        disable join of overlapping reads [default: %d]\n", p_join);
    fprintf(stderr, "    -c        disable chop of overlapping reads [default: %d]\n\n", p_chop);
    fprintf(stderr, "The default settings:\n");
    fprintf(stderr, " - chop between  6- 8bp with 0 mismatches\n");
    fprintf(stderr, " - chop between  9-13bp with 1 mismatch\n");
    fprintf(stderr, " - join between 14-16bp with 1 mismatch\n");
    fprintf(stderr, " - join between 17-24bp with 2 mismatches\n");
    fprintf(stderr, " - join between 25-  bp with 3+ mismatches\n");
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
    int h[3] = {0}; // 
    int last0, last1; // index of last base
    char *s0, *s1; // sequence
    char *q1; // quality scores 
    char *rc1, *rq1; // reverse-complement and reverse-qualities of s1
    rc1 = calloc(MAX_READ_LENGTH * sizeof(char), sizeof(char));
    rq1 = calloc(MAX_READ_LENGTH * sizeof(char), sizeof(char));
    FILE *fw[3];
    gzFile  fp[2];
    kseq_t *ks[2];
    kseq_t *r;
    FILE *f;
    // file io
    fw[0] = fopen(fo[0], "wb");
    fw[1] = fopen(fo[1], "wb");
    fw[2] = fopen(fo[2], "wb");
    fp[0] = gzopen(argv[optind], "rb");
    fp[1] = gzopen(argv[optind+1], "rb");
    ks[0] = kseq_init(fp[0]);
    ks[1] = kseq_init(fp[1]);
    int operation; // PASS, CHOP, or JOIN
    int min_off = p_min_off_chop; // guaranteed lower than p_min_off_join
    double max_err = p_max_err_chop; // guaranteed higher than p_max_err_join
    while ((kseq_read(ks[0]) >= 0) && (kseq_read(ks[1]) >= 0)) {
      bt = 0;
      operation = PASS;
      last1 = ks[1]->seq.l-1;
      s1 = ks[1]->seq.s;
      for (i=0; i<=last1; ++i) {
        // reverse complement sequence
        switch(s1[i]) {
        case 'A': rc1[last1 - i] = 'T'; break;
        case 'T': rc1[last1 - i] = 'A'; break;
        case 'G': rc1[last1 - i] = 'C'; break;
        case 'C': rc1[last1 - i] = 'G'; break;
        default:  rc1[last1 - i] = 'N'; break;
        }
      }
      s1 = rc1;
      last0 = ks[0]->seq.l-1;
      s0 = ks[0]->seq.s;
      //
      int i = (last0 < last1 ? last0 : last1) + 1;
      while (i > min_off){
        i--;
        h[0] = HammingDistance(&s0[last0 - (i - 0)], (i-0)+1, (int)((i+1) * max_err), s1);
        if (p_join && (i >= p_min_off_join) && (h[0] <= ((i+1) * p_max_err_join))) {
          max_err = p_max_err_join;
          operation = JOIN;
          break;
        }
        if (p_chop && (i >= p_min_off_chop) && (h[0] <= ((i+1) * p_max_err_chop))) {
          max_err = p_max_err_chop;
          operation = CHOP;
          break;
        }
      }
      /* for (i=offset; i>=min_off; i--) { */
      /* } */
      //
      if (operation != PASS) {
        // calculate two backtracks
        h[1] = HammingDistance(&s0[last0 - (i-1)], (i-1) + 1, (int)(((i-1)+1) * max_err), s1);
        h[2] = HammingDistance(&s0[last0 - (i-2)], (i-2) + 1, (int)(((i-2)+1) * max_err), s1);
        // <= because we favor lower backtrack
        bt = (h[0] <= h[1]) ? (h[0] <= h[2] ?   0  :   2 ) : (h[1] <= h[2] ?   1 :    2 );
        i-=bt;
      }
      // 
      opcounts[operation]++;
      if (verbose > 0) {
        int hd=0; // hamming distance
        if (operation != PASS) {
          hd = (h[0] <= h[1]) ? (h[0] <= h[2] ? h[0] : h[2]) : (h[1] <= h[2] ? h[1] : h[2]);
        } else {
          hd = h[0];
        }
        printf("operation:%d ", operation);
        printf("overlap:%d dinstance:>=%d backtrack:%d \n%s\n", i+1, hd, bt, ks[0]->seq.s);
        for (int z=0; z < (last0-i); z++) {
          printf(" ");
        }
        printf("%s\n\n", rc1);
         verbose--;
      }
      
      if (operation != PASS) {
        ks[0]->seq.s[last0 - i] = '\0';
        ks[0]->qual.s[last0 - i] = '\0';
      }

      if (operation == JOIN) {
        // reverse quality scores
        q1 = ks[1]->qual.s;
        for (i=0; i<=last1; ++i) {
          rq1[last1 - i] = q1[i];
        }
        f = fw[0];
        fputc('@', f);
        fputs(ks[0]->name.s, f);
        fputc('\n', f);
        fputs(ks[0]->seq.s, f);
        fputs(rc1, f);
        fputc('\n', f);
        fputs("+", f);
        fputc('\n', f);
        fputs(ks[0]->qual.s, f);
        fputs(rq1, f);
        fputc('\n', f);
      } else {
        //
        f = fw[1];
        r = ks[0];
        fputc('@', f);
        fputs(r->name.s, f);
        fputc('\n', f);
        fputs(r->seq.s, f);
        fputc('\n', f);
        fputs("+", f);
        fputc('\n', f);
        fputs(r->qual.s, f);
        fputc('\n', f);
        f = fw[2];
        r = ks[1];
        fputc('@', f);
        fputs(r->name.s, f);
        fputc('\n', f);
        fputs(r->seq.s, f);
        fputc('\n', f);
        fputs("+", f);
        fputc('\n', f);
        fputs(r->qual.s, f);
        fputc('\n', f);
      }
    }
    // cleanup
    free(rc1);
    free(rq1);
    kseq_destroy(ks[0]);
    kseq_destroy(ks[1]);
    gzclose(fp[0]);
    gzclose(fp[1]);
    fclose(fw[0]);
    fclose(fw[1]);
    fclose(fw[2]);
  }
  // summary
  printf("JOIN'ed:\t%" PRIu64 "\n", opcounts[JOIN]);
  printf("CHOP'ed:\t%" PRIu64 "\n", opcounts[CHOP]);
  printf("PASS'ed:\t%" PRIu64 "\n", opcounts[PASS]);
  // cleanup
  free(fo[0]);
  free(fo[1]);
  free(fo[2]);
  return status;
}
