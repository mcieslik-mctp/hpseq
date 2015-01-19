#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"
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

typedef struct {
        int type, len;
        char *seq;
        uint64_t cnt;
} ta_adap_t;

int main(int argc, char *argv[]) {
  int n_adaps, m_adaps;
  int from_stdin;
  int c, j, i;

  double max_err = 0.1;
  int min_cut = 3;
  int min_len = 25;
  int trim = 1;
  int wipe = 1;
  
  ta_adap_t *adaps;
  kseq_t *ks;
  gzFile fp;

  n_adaps = m_adaps = 0; adaps = 0;
  while ((c = getopt(argc, argv, "5:3:e:c:l:t:w:")) >= 0) {
    if (c == '5' || c == '3') {
      ta_adap_t *p;
      if (m_adaps == n_adaps) {
	m_adaps = m_adaps? m_adaps<<1 : 4;
	adaps = realloc(adaps, m_adaps * sizeof(ta_adap_t));
      }
      p = &adaps[n_adaps++];
      p->seq = (char*)strdup(optarg);
      p->type = c - '0';
    }
    else if (c == 'e') max_err = atof(optarg);
    else if (c == 'c') min_cut = atoi(optarg);
    else if (c == 'l') min_len = atoi(optarg);
    else if (c == 't') trim = atoi(optarg);
    else if (c == 'w') wipe = atoi(optarg);
  }
  
  // preset
  if (n_adaps == 0) {
    m_adaps = n_adaps = 2;
    adaps = malloc(m_adaps * sizeof(ta_adap_t));
    adaps[0].type = 3;
    adaps[0].seq = (char*)strdup("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
    adaps[1].type = 3;
    adaps[1].seq = (char*)strdup("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT");
  }
  // update adapter info
  for (j = 0; j < n_adaps; ++j) {
    ta_adap_t *p = &adaps[j];
    p->len = strlen((char*)p->seq);
    p->cnt = 0;
  }

  from_stdin = !isatty(fileno(stdin));
  if (optind == argc && !from_stdin) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   thecut [options] <in.fq>\n\n");
    fprintf(stderr, "         -3 STR     3'-end adapter\n");
    fprintf(stderr, "         -e INT     max adapter alignment error [%.3f]\n", max_err);
    fprintf(stderr, "         -l INT     min length of remaining sequence [%d]\n", min_len);
    fprintf(stderr, "         -c INT     min length of cut adapter [%d]\n", min_cut);
    fprintf(stderr, "         -t INT     if 0 no triming is done [%d]\n", trim);
    fprintf(stderr, "         -w INT     if 0 no wiping is done [%d]\n", wipe);
    fprintf(stderr, "\n");
    return 1; // FIXME: memory leak
  }

  fp = optind < argc && strcmp(argv[optind], "-")? gzopen(argv[optind], "rb") : gzdopen(fileno(stdin), "rb");
  ks = kseq_init(fp);

  uint64_t counts[200] = {0};
  while (kseq_read(ks) >= 0) {
    int end=ks->seq.l;
    int l;
    int hd;
    int e;
    int mm;
    for (j = 0; j < n_adaps; ++j) {
      ta_adap_t *p = &adaps[j];
      if (p->type == 3) {
	for (i=ks->seq.l; i>=min_cut; --i) {
	  l = i > p->len ? p->len : i;
	  e = ks->seq.l - i;
	  mm = (int)(l * max_err);
	  hd = HammingDistance(p->seq, l, mm, &ks->seq.s[e]);
	  if (hd <= mm) {
	    p->cnt++;
	    if (e < end) {
	      end=e;
	      break;
	    }
	  }
	}
      }
    }
    counts[ks->seq.l - end]++;
    if (wipe) {
      for (i = end; i < ks->seq.l; ++i) {
	ks->seq.s[i] = 'N';
      }
    }

    putchar(ks->qual.l? '@' : '>');
    puts(ks->name.s);
    if (!trim | (end < min_len)) {
      puts(ks->seq.s);
    } else {
      fwrite(ks->seq.s, 1, end, stdout); putchar('\n');
    }
    if (ks->qual.l) {
      puts("+");
      if (!trim | (end < min_len)) {
  	puts(ks->qual.s);
      } else {
  	fwrite(ks->qual.s, 1, end, stdout); putchar('\n');
      }
    }
  }
  
  kseq_destroy(ks);
  gzclose(fp);
  
  for (j = 0; j < n_adaps; ++j) {
    ta_adap_t *p = &adaps[j];
    fprintf(stderr, "adapter:\t%s\t%ld\n", p->seq, (long)p->cnt);
    free(p->seq);
  }

  for (j = 1; j < 20; j++) {
    fprintf(stderr, "length:\t%d\t%d\n", j, (int)counts[j]);
  }
  
  free(adaps);
  return 0;
}
