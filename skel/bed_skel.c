#include <stdlib.h>
#include <stdio.h>
#include <hts.h>
#include <sam.h>
#include <string.h>

int main(int argc, char *argv[]) {
  int exit_code = 0;
  // io names and handles
  char *out_bed_fn;
  out_bed_fn = malloc(strlen(argv[1])+1);
  strcpy(out_bed_fn, argv[1]);
  strcpy(strrchr(out_bed_fn, '.'), ".bed");
  htsFile *inp;
  FILE    *out_bed;

  int ret;
  bam1_t *br;
  bam_hdr_t *hdr;
  const htsFormat *fmt;

  //
  inp = hts_open(argv[1], "r");
  fmt = hts_get_format(inp);
  hdr = sam_hdr_read(inp);
  out_bed = fopen(out_bed_fn, "wb");
  fprintf(out_bed, "# format: %d\n", fmt->format);
  br = bam_init1();
  int i=0;
  while ((ret = sam_read1(inp, hdr, br)) >= 0) {
    fprintf(out_bed, "%s\t%d\t%d\t%c\n",
	    hdr->target_name[br->core.tid], br->core.pos,
  	    bam_endpos(br), bam_is_rev(br) ? '+' : '-');
    i++;
  }
  bam_destroy1(br);
  bam_hdr_destroy(hdr);
  fclose(out_bed);
  free(out_bed_fn);

  ret = sam_close(inp);
  if (ret < 0) {
    fprintf(stderr, "Error closing input.\n");
    exit_code = 1;
  }
  
  return exit_code;
}
