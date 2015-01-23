#include <stdlib.h>
#include <stdio.h>
#include <hts.h>

int main() {

  htsFile* hfi = hts_open("test/pe_alig.bam", "r");
  const htsFormat* z;
  z = hts_get_format(hfi);
  printf("%d\n", z->format);
  
  return 0;
}
