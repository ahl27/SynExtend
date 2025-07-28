#ifndef EXOLABEL_H
#define EXOLABEL_H

#ifdef COMPILING_SYNEXTEND_VIA_R
  #include "../SEutils.h"
#else
  #include "FallbackDefines.h"
#endif

#include "PrefixTrie.h"
#include "LoserTree.h"
#include "FileHandlers.h"
#include <time.h>

// includes for file truncation
#ifdef HAVE_UNISTD_H
  #include <unistd.h>
#else
  #ifdef WIN32
    #include <io.h>
  #endif
#endif

/***********/
/* Defines */
/***********/
#define h_uint uint64_t
#define uint uint_fast32_t
#define strlen_uint uint_least16_t
#define l_uint uint64_t
#define lu_fprint PRIu64
#define w_float float
#define aq_int int16_t // size of "seen" counter in ArrayQueue

/*
 * common limits are defined in limits.h
 *  PAGE_SIZE: size in bytes of a page
 *  PATH_MAX: max size in bytes of a path name
 *  I'm not really using these for anything yet, but maybe eventually
 */
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#ifndef PAGE_SIZE
#define PAGE_SIZE 4096
#endif

// max size of a vertex name (char array will have one extra space for terminator)
#define MAX_NODE_NAME_SIZE 255

// holds pointers char* of size MAX_NODE_NAME_SIZE, 4096 is 1MB
#define NODE_NAME_CACHE_SIZE 40960

// number of entries, so total consumption often multiplied by 8
#define FILE_READ_CACHE_SIZE 8192*4


/*************/
/* Constants */
/*************/

// Constants for weight compression, see compressEdgeWeights for details
static const int BITS_FOR_WEIGHT = 16;
static const int BITS_FOR_EXP = 4;

static const double MIN_DOUBLE = -1 * DBL_MIN;
static const int L_SIZE = sizeof(l_uint);
static const int W_SIZE = sizeof(w_float);

/*
 * Some comments on external sorting performance:
 *  - FILE_READ_CACHE_SIZE determines the initial sort size and the buffer size
 *  - More buffers means more RAM consumption, but fewer passes through the file
 *  - Block interchanges happen when some not-yet-processed data will be
 *    overwritten by an output dump. These are the slowest operation. Could maybe
 *    be sped up by using a fixed size temporary file as a larger buffer...?
 *    More intelligent block interchanges would also speed up a lot.
 *  - Larger input buffers means more sequential access and fewer block
 *    interchange operations
 *  - More buffers means more block interchange operations
 *  - Mergesort bin space is multiplied by sizeof(edge) (>=16 bytes)
 */
static const int MAX_BINS_FOR_MERGE = 64; // will round up to next highest power of 2
static const int MERGE_INPUT_SIZE = FILE_READ_CACHE_SIZE;
static const int MERGE_OUTPUT_SIZE = 16*FILE_READ_CACHE_SIZE;

// Erik says this isn't really useful, so disabling
// set to a positive value for a cutoff to re-add values to the queue
static const int MAX_EDGES_EXACT = -1;

// Progress printing values (also controls R_checkUserInterrupt() checks)
static const int PRINT_COUNTER_MOD = 811*13;
static const int PROGRESS_COUNTER_MOD = 6043;


/**********************/
/* Exported functions */
/**********************/
void cleanup_ondisklp_global_values();

void cluster_file(const char* weights_fname,
                  const char* neighbor_fname,
                  const l_uint num_v, const int max_iterations, const int v,
                  const float self_loop_weight,
                  const double atten_pow);

#endif