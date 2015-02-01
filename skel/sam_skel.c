#include <stdlib.h>
#include <stdio.h>
#include <hts.h>
#include <sam.h>
#include <string.h>


int main(int argc, char *argv[]) {
  int exit_code = 0;
  // io names and handles
  char *out_sam_fn;
  out_sam_fn = malloc(strlen(argv[1])+1);
  strcpy(out_sam_fn, argv[1]);
  strcpy(strrchr(out_sam_fn, '.'), ".sam");
  htsFile *inp;
  htsFile *out_sam;
  // 
  int ret;
  bam1_t *br;
  bam_hdr_t *hdr;

  const htsFormat *fmt;

  inp = hts_open(argv[1], "r");
  hdr = sam_hdr_read(inp);
  out_sam = hts_open(out_sam_fn, "w");
  sam_hdr_write(out_sam, hdr);

  br = bam_init1();
  while ((ret = sam_read1(inp, hdr, br)) >= 0) {
    // do something to br
    if (sam_write1(out_sam, hdr, br) < 0) {
      fprintf(stderr, "Error writing output.\n");
      exit_code = 1;
      break;
    }
  }
  ret = sam_close(out_sam);
  if (ret < 0) {
    fprintf(stderr, "Error closing output.\n");
    exit_code = 1;
  }
  ret = sam_close(inp);
  if (ret < 0) {
    fprintf(stderr, "Error closing input.\n");
    exit_code = 1;
  }
  bam_destroy1(br);
  bam_hdr_destroy(hdr);
  free(out_sam_fn);
  
  return exit_code;
}
