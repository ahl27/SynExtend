#ifndef COMPILING_SYNEXTEND_VIA_R
/*
 * Fallback definitions to be used if not compiling within R
 * e.g., if using as a standalone
 * COMPILING_SYNEXTEND_VIA_R is set in Makevars in the R installation
 */

  #ifndef FALLBACK_NONR_DEFINES
  #define FALLBACK_NONR_DEFINES

  #include <stdio.h>
  #include <stdlib.h>
  #include <stdarg.h>
  #include <math.h>
  #include <float.h>
  #include <time.h>
  #include <inttypes.h>

  #define safe_malloc malloc
  #define safe_calloc calloc
  #define safe_realloc realloc


  // Standalone fallback for Rprintf
  static inline void Rprintf(const char *fmt, ...) {
      va_list args;
      va_start(args, fmt);
      vprintf(fmt, args);
      va_end(args);
  }

  static inline void error(const char *fmt, ...) {
      va_list args;
      va_start(args, fmt);
      vfprintf(stderr, fmt, args);
      va_end(args);
      fprintf(stderr, "\n");
      exit(EXIT_FAILURE);
  }

  static inline uint32_t get_msb32(uint32_t v){
    // this is taken from https://stackoverflow.com/a/31718095
    static const int debruijn_bitmult[32] =
    {
        0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return debruijn_bitmult[(uint32_t)(v * 0x07C4ACDDU) >> 27];
  }

  static inline void R_CheckUserInterrupt(void) { return; }
  static inline void GetRNGstate(void) { srand(time(NULL)); }
  static inline void PutRNGstate(void) { return; }
  static inline double unif_rand(void) { return ((double)rand()) / RAND_MAX; }

  #endif
#endif