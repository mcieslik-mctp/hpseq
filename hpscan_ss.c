#include <stdio.h>
#include <emmintrin.h>
#include <unistd.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>

static void setbit(void *v, int p) { ((uint32_t*)v)[p >> 5] |= 1 << (p & 31); }

char *bndm128(char *target, int tgtlen, char *pattern, int patlen)
{
  assert(patlen <= 128);
  uint8_t     *tgt = (uint8_t*)target, *pat = (uint8_t*)pattern;
  int         i, j;
  __m128i     zero = {0}, maskv[256], carry64 = (__v2di){0,1};
  uint8_t     used[256] = {0};
  for (i = 0; i < patlen; ++i) {
    if (!used[j = pat[i]])
      used[j] = 1, maskv[j] = zero;
    setbit(&maskv[j], patlen - 1 - i);
  }

  for (i = 0; i <= tgtlen - patlen; i += j) {
    j = patlen;
    if (!used[tgt[patlen - 1 + i]]) continue;
    __m128i     mask = maskv[tgt[patlen - 1 + i]];
    while (0xFFFF != _mm_movemask_epi8(_mm_cmpeq_epi8(zero, mask))) {
      if (!--j)
	return target + i;
      if (!used[tgt[i + j - 1]])
	break;
      mask = _mm_or_si128( _mm_slli_epi64(mask, 1)
			   , _mm_srli_epi64( _mm_slli_si128(mask, 8), 63 )
			   );
    }
  }
  return 0;
}

int main(int argc, char *argv[]) {
}
