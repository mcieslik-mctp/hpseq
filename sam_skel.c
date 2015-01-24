#include <stdlib.h>
#include <stdio.h>
#include <hts.h>
#include <sam.h>


int main(int argc, char *argv[]) {
  int ret;
  int exit_code = 0;
  const htsFormat* fmt;
  bam_hdr_t *hdr;
  htsFile* inp;
  htsFile* out;
  inp = hts_open(argv[1], "r");
  out = hts_open("-", "w");

  fmt = hts_get_format(inp);
  printf("%d\n", fmt->format);
  
  hdr = sam_hdr_read(inp);
  sam_hdr_write(out, hdr);

  //
  bam1_t *b;
  b = bam_init1();
  while ((ret = sam_read1(inp, hdr, b)) >= 0) {
    if (sam_write1(out, hdr, b) < 0) {
      fprintf(stderr, "Error writing output.\n");
      exit_code = 1;
      break;
    }
  }
  // clean-up
  bam_destroy1(b);
  bam_hdr_destroy(hdr);
  
  ret = sam_close(out);
  if (ret < 0) {
    fprintf(stderr, "Error closing output.\n");
    exit_code = 1;
  }
  ret = sam_close(inp);
  if (ret < 0) {
    fprintf(stderr, "Error closing input.\n");
    exit_code = 1;
  }
  
  return exit_code;
}
