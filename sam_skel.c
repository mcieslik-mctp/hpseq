#include <stdlib.h>
#include <stdio.h>
#include <hts.h>
#include <sam.h>
#include <string.h>


int main(int argc, char *argv[]) {
  int ret;
  int exit_code = 0;
  bam_hdr_t *hdr;
  htsFile *inp;
  htsFile *out_sam;
  FILE *out_bed;
  const htsFormat *fmt;
  
  inp = hts_open(argv[1], "r");
  hdr = sam_hdr_read(inp);
  fmt = hts_get_format(inp);
  fprintf(stderr, "format: %d\n", fmt->format);

  
  // change extension of output file
  strcpy(strrchr(argv[1], '.'), ".sam");
  out_sam = hts_open(argv[1], "w");
  sam_hdr_write(out_sam, hdr);

  // change extension of output file
  strcpy(strrchr(argv[1], '.'), ".bed");
  // 
  out_bed = fopen(argv[1], "wb");
  bam1_t *br;
  br = bam_init1();
  while ((ret = sam_read1(inp, hdr, br)) >= 0) {
    fprintf(out_bed, "%s\t%d\t%d\t%c\n", hdr->target_name[br->core.tid], br->core.pos,
	    bam_endpos(br), bam_is_rev(br) ? '+' : '-');
  }
  bam_destroy1(br);
  fclose(out_bed);

  br = bam_init1();
  rewind(inp);
  while ((ret = sam_read1(inp, hdr, br)) >= 0) {
    if (sam_write1(out_sam, hdr, br) < 0) {
      fprintf(stderr, "Error writing output.\n");
      exit_code = 1;
      break;
    }
  }
  // clean-up
  bam_destroy1(br);
  bam_hdr_destroy(hdr);
  
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
  
  return exit_code;
}
