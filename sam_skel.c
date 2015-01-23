#include <stdlib.h>
#include <stdio.h>
#include <hts.h>
#include <sam.h>


int main(int argc, char *argv[]) {
  const htsFormat* fmt;
  bam_hdr_t *hdr;
  htsFile* inp;
  /* htsFile* out; */
  inp = hts_open(argv[1], "r");
  /* out = hts_open("-", "w"); */

  fmt = hts_get_format(inp);
  printf("%d\n", fmt->format);
   
  /*  sam_index_load */
  /* sam_hdr_write(out, h); */
  return 0;
}
