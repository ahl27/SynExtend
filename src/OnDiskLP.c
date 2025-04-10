/*
 * Out of memory clustering with fast label propagation
 * Author: Aidan Lakshman
 *
 *
 * Vertices are stored in a trie. Trie leaves have the following attributes:
 *  -      count: Initially stores the degree of the node, later its cluster
 *  -      index: Stores the index of the node. May be able to remove later.
 *  - edge_start: Stores the offset to relevant edges in the below files
 *  -       dist: Stores the distance used for attenuation
 *
 *  This creates a csr-compressed graph structure. The index is stored in the
 *  `edge_start` values of trie leaves, which is a total of n+1 uint64_t values.
 *  The k'th value denotes the start position in the below files for vertex k.
 *  For example, if the first two values are 0 100, then the incoming edges to
 *  the first vertex are entries 0-99.
 *
 * A total of 3 files are created:
 *  - neighbors: source node for each incoming edge (indexed by trie).
 *  -   weights: weight for each incoming edge (indexed by trie).
 *  -   outfile: .tsv file returned to R, contains two tab-separated columns
 *               (vertex name, cluster)
 *
 * Additional Notes:
 *  - sizeof(char) is guaranteed to be 1 (see https://en.wikipedia.org/wiki/Sizeof).
 *    1 is used instead to simplify code somewhat.
 *  - fwrite is thread safe but gives no performance benefit to parallelize
 *    (https://stackoverflow.com/questions/26565498/multiple-threads-writing-on-same-file)
 *
 * TODOs:
 *  - At some point it's probably worth refactoring all file accesses into some
 *    kind of struct w/ accessors
 *    -> something like a virtual array object that's actually r/w to disk,
 *        could be useful in future
 *    -> this implementation should use mmap (and Windows equivalent when
 *        necessary) to improve random r/w
 *    -> this is a very far off wishlist item, not necessary
 *
 * Optimization Ideas:
 *  - Can we speed up reading edges more?
 *    - add a cache of "recently used vertices" as [c1,c2,...,\0,ptr]
 *    - is it possible to order the nodes by what is most often referenced?
 *      ...seems like more trouble than its worth
 *    - reduce number of block interchanges in the in-place mergesort
 *  - Can we optimize RAM consumption?
 *    - might be worthwhile to compress the trie structure *after* its read in
 *    - can probably remove the `index` attribute of trie leaves
 */

#include "SEutils.h"
#include "SynExtend.h"
#include "PrefixTrie.h"
#include "LoserTree.h"
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

static const double MIN_DOUBLE = -1 * DBL_MIN;
static const int L_SIZE = sizeof(l_uint);
static const int W_SIZE = sizeof(w_float);
static const char CONSENSUS_CSRCOPY1[] = "tmpcsr1";
static const char CONSENSUS_CLUSTER[] = "tmpclust";

/*
 * Constants for weight compression, see compressEdgeWeights for details
 */
static const int BITS_FOR_WEIGHT = 16;
static const int BITS_FOR_EXP = 4;

// don't touch, values auto-determined from above constants
static const l_uint MAX_NUM_NODES = (1ULL << (64-(BITS_FOR_WEIGHT+BITS_FOR_EXP))) - 1;
static const l_uint MAX_POSSIBLE_WEIGHT = (1ULL << BITS_FOR_WEIGHT) - 1;
static const l_uint MAX_POSSIBLE_EXPONENT = (1ULL << BITS_FOR_EXP) - 1;

enum VERBOSITY {
  VERBOSE_NONE = 0,   // no output
  VERBOSE_BASIC = 1,  // non-interactive output, better for files
  VERBOSE_ALL = 2     // interactive output, lots of carriage returns
};

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
/* Struct Definitions */
/**********************/
typedef struct {
  l_uint *queue;
  aq_int *seen; // can make this bigger if necessary
  aq_int max_seen;
  l_uint start;
  l_uint end;
  l_uint size;
  l_uint length;
} ArrayQueue;

typedef struct {
  l_uint v;
  l_uint w;
} edge;

/*********************************/
/* Global Vars to Cleanup at End */
/*********************************/
static leaf **GLOBAL_leaf = NULL;
static leaf **GLOBAL_all_leaves = NULL;
static char **GLOBAL_filenames = NULL; // holds file NAMES
static int GLOBAL_nfiles = 0;
static l_uint CLUST_MAP_CTR = 0;
static ArrayQueue *GLOBAL_queue = NULL;
static prefix *GLOBAL_trie = NULL;
static FILE **GLOBAL_ftracker = NULL; // holds file POINTERS
static int GLOBAL_num_files = 0;
static int GLOBAL_nbuffers = 0;
static void **GLOBAL_mergebuffers = NULL;
static LoserTree *GLOBAL_mergetree = NULL;
static l_uint GLOBAL_verts_changed = 0;
static edge* GLOBAL_readedges = NULL;
static l_uint GLOBAL_cachectr = 0;

/***************************/
/* Struct Helper Functions */
/***************************/

static void free_array_queue(ArrayQueue *q){
  free(q->queue);
  free(q->seen);
  free(q);
}

static void array_queue_insert(ArrayQueue *q, l_uint value){
  // if the value is already in queue, just return
  // this is denoted by having a positive value
  // (0 denotes "we've seen this vertex enough")
  int nseen = q->seen[value];
  if(nseen >= 0) return;

  // else, add to queue:
  // find next position and increment end
  l_uint pos = q->end;
  q->end = (pos+1) % q->size;
  q->queue[pos] = value;

  // flip sign of nseen to positive to add to queue
  // increment nseen by one, set to zero if we're done
  nseen = (-1*nseen) + 1;
  if(nseen == q->max_seen) nseen = 0;
  q->seen[value] = nseen;
  q->length++;
  return;
}

static l_uint array_queue_pop(ArrayQueue *q){
  if(!q->length) error("Attempted to pop from queue with no elements.");
  l_uint val = q->queue[q->start];
  q->seen[val] *= -1; // set to negative to denote "removed from queue"
  q->start = (q->start + 1) % q->size;
  q->length--;
  return val;
}

static ArrayQueue* alloc_array_queue(l_uint size, int max_seen){
  if(!size) error("Attempted to initialize queue of size 0.");
  ArrayQueue *q = safe_malloc(sizeof(ArrayQueue));
  q->queue = safe_malloc(L_SIZE * size);
  q->seen = safe_malloc(sizeof(aq_int) * size);
  q->max_seen = max_seen;
  q->start = 0;
  q->end = 0;
  q->length = 0;
  q->size = size;
  return q;
}

static ArrayQueue* init_array_queue(l_uint size, int max_seen){
  ArrayQueue *q = alloc_array_queue(size+1, max_seen);

  // initialize with a random permutation
  GetRNGstate();
  l_uint j;
  for(l_uint i=0; i<size; i++){
    j = (l_uint) trunc((i+1) * (unif_rand()));
    q->queue[i] = i;
    if(j < i){
      q->queue[i] = q->queue[j];
      q->queue[j] = i;
    }
    q->seen[i] = 1; // set each vertex to have been seen once
  }
  PutRNGstate();

  // guard case: we have an extra slot, mark it as a never visit
  q->seen[size] = 0;

  // now start and end are both 0
  // this is fine, since end is the next location to insert to and queue is full
  q->end = size;
  q->length = size;
  return q;
}

/********************/
/* Output Functions */
/********************/

static void report_time(time_t time1, time_t time2, const char* prefix){
  double elapsed_time = difftime(time2, time1);
  int mins, hours, days, secs;
  secs = (int)fmod(elapsed_time, 60);

  l_uint elapsed_secs = (l_uint)(elapsed_time - secs);
  days = elapsed_secs / 86400;
  elapsed_secs %= 86400;

  hours = elapsed_secs / 3600;
  elapsed_secs %= 3600;

  mins = elapsed_secs / 60;

  Rprintf("%sTime difference of ", prefix);
  if(days) Rprintf("%d day%s, ", days, days==1 ? "" : "s");
  if(hours) Rprintf("%d hr%s, ", hours, hours==1 ? "" : "s");
  if(mins) Rprintf("%d min%s, ", mins, mins==1 ? "" : "s");
  Rprintf("%d sec%s\n", secs, secs==1 ? "" : "s");
  return;
}

static void report_filesize(l_uint bytes){
  // this will assume the file pointer is already at the end
  double nbytes = (double)bytes;
  char units[] = {' ', 'K', 'M', 'G', 'T', 'P'};
  char unit_out[] = " B";
  int cur_pos = 0;
  while(nbytes > 1000 && cur_pos < 6){
    nbytes /= 1000;
    cur_pos++;
  }
  unit_out[0] = units[cur_pos];
  Rprintf("%5.1f%s", nbytes, unit_out);
  return;
}

static void print_graph_stats(l_uint num_v, l_uint num_e){
  l_uint vals[] = {num_v, num_e};
  l_uint divisor;
  int to_print, first_val, is_really_big;
  for(int i=0; i<2; i++){
    is_really_big = vals[i] >= (1000000000ULL * (vals[i] ? 1 : 1000));
    if(i == 0)
      Rprintf("Total Vertices: ");
    else
      Rprintf("Total Edges: ");
    divisor = 1;
    first_val = 1;
    while(divisor*1000 <= vals[i]) divisor *= 1000;
    while(divisor > 1){
      to_print = vals[i] / divisor;
      if(first_val){
        Rprintf("%d,", to_print);
        first_val = 0;
      } else Rprintf("%03d,", to_print);
      vals[i] %= divisor;
      divisor /= 1000;
    }
    to_print = vals[i];
    if(first_val){
        Rprintf("%d%s\n", to_print, is_really_big ? " (WOW!)" : "");
        first_val = 0;
    } else Rprintf("%03d%s\n", to_print, is_really_big ? " (WOW!)" : "");
  }
}

/************************/
/* Arithmetic Functions */
/************************/

static void kahan_accu(double *cur_sum, double *cur_err, double new_val){
  // this is an accumulator that uses kahan summation to improve floating point accuracy
  double sum, y, z;
  sum = *cur_sum;

  // add in the previous error
  y = *cur_err + new_val;

  // now we do Fast2Sum(a=sum, b=y)
  *cur_sum = sum + y;
  z = *cur_sum - sum;
  *cur_err = y - z;

  return;
}

static inline float sigmoid_transform(const float w, const double slope){
  // should probably expose these at some point
  const float scale = 0.5;
  const float cutoff = 0.1;
  float r = 1 / (1+exp(-1*slope*(w-scale)));
  return r > cutoff ? r : 0;
}

static edge compressEdgeValues(l_uint v1, l_uint v2, double weight){
  /*
   * Using a compressed representation for weight
   * BITS_FOR_WEIGHT determines the number of bits for the value
   * BITS_FOR_EXP the number of bits for the exponent.
   * The compression follows the following algorithm:
   *  1. Transform each weight w with w' = log2(1+w).
   *  2. Record the position of the leading bit of the integer part, M
   *      Note M = min(M, 1ULL << BITS_FOR_EXP)
   *  3. Right-shift w' by (BITS_FOR_WEIGHT - M)
   *      (multiply by (1ULL << (BITS_FOR_WEIGHT - M)))
   *  4. Store the value in the lower [BITS_FOR_WEIGHT] bits
   *
   * Note that weights are guaranteed to be positive.
   * Performance:
   * - Maximum error of 0.00004 for seq(0,1,0.001)
   * - Maximum error of around 10 for >60,000 (~0.01%)
   * - Basically unbounded range of values
   * - Absolute error goes up a lot as we increase, but typically not more than 0.05%
   */
  if(v2 > MAX_NUM_NODES){
    cleanup_ondisklp_global_values();
    error("ExoLabel can only support up to %llu nodes.", MAX_NUM_NODES);
  }
  edge e;
  e.v = v1;

  // left shift the index out of the range used for the weight
  l_uint v2_comp = v2;
  v2_comp <<= (BITS_FOR_WEIGHT + BITS_FOR_EXP);

  // 1. transform the weight into log2(w+1)
  weight = log2(weight+1);
  if(weight > MAX_POSSIBLE_WEIGHT) weight = MAX_POSSIBLE_WEIGHT;

  // 2. get the position of the leading bit
  uint32_t to_shift = MAX_POSSIBLE_EXPONENT - get_msb32((uint32_t)floor(weight));
  if(to_shift < 0) to_shift = 0;

  // 3. right-shift w' by (BITS_FOR_WEIGHT - M)
  weight *= (1ULL << to_shift);

  // mask the values to make sure they're in the right position
  l_uint w_comp = ((l_uint)floor(weight)) & MAX_POSSIBLE_WEIGHT;
  to_shift = (to_shift & MAX_POSSIBLE_EXPONENT) << BITS_FOR_WEIGHT;

  e.w = v2_comp | to_shift | w_comp;
  return e;
}

static void decompressEdgeValue(l_uint compressed, l_uint *v2, w_float *w){
  /*
   * See compressEdgeValues for description
   * 64-bit unsigned integer
   * first 48 bits are the index
   * last 16 bits are the compressed weight (4.12 format)
   */

  // get the weight in compressed form
  double w_comp = (double) ((l_uint)(compressed & MAX_POSSIBLE_WEIGHT));

  // get the exponent that we shifted by
  uint32_t bits_shifted = (compressed >> BITS_FOR_WEIGHT) & MAX_POSSIBLE_EXPONENT;

  w_comp /= (1ULL << bits_shifted);
  w_comp = pow(2, w_comp) - 1;

  *w = (w_float)w_comp;
  *v2 = compressed >> (BITS_FOR_WEIGHT + BITS_FOR_EXP);

  return;
}

static void safe_filepath_cat(const char *dir, const char *f, char *fname, size_t fnamesize){
  // fname should be preallocated
  char directory_separator;
#ifdef WIN32
  directory_separator = '\\';
#else
  directory_separator = '/';
#endif
  memset(fname, 0, fnamesize);
  // length is len(dir) + len(f) + 1 (separator) + 1 (terminator)
  snprintf(fname, strlen(dir)+strlen(f)+2, "%s%c%s", dir, directory_separator, f);
  return;
}

static FILE* safe_fopen(const char *filename, const char *mode){
  FILE *f = fopen(filename, mode);
  GLOBAL_ftracker[GLOBAL_num_files++] = f;
  return f;
}

static void fclose_tracked(int nfiles){
  if(!nfiles) return;
  if(nfiles > GLOBAL_num_files) error("attempted to close more files than were open!");
  // going to use a PROTECT-like strategy
  // just close the top N files
  FILE *f;
  for(int i=0; i<nfiles; i++){
    // check to make sure we dont fclose(NULL)
    f = GLOBAL_ftracker[--GLOBAL_num_files];
    if(f) fclose(f);
  }
  return;
}

static char* create_filename(const char* dir, const char* to_append){
  size_t nchar1 = strlen(to_append);
  size_t nchar2 = strlen(dir) + nchar1 + 3;
  char *newstr = safe_malloc(nchar2);
  safe_filepath_cat(dir, to_append, newstr, nchar1);
  // save the result somewhere for later
  GLOBAL_filenames[GLOBAL_nfiles++] = newstr;
  return newstr;
}

void cleanup_ondisklp_global_values(){
  /*
   * This is called with on.exit
   * note that it could be called even if the main function never runs!
   * be cautious with freeing pointers that may be NULL
   */

  // cleanup open file handles
  fclose_tracked(GLOBAL_num_files);
  if(GLOBAL_ftracker){
    free(GLOBAL_ftracker);
    GLOBAL_ftracker = NULL;
  }

  // free array queue
  if(GLOBAL_queue){
    free_array_queue(GLOBAL_queue);
    GLOBAL_queue = NULL;
  }

  // delete any files made during runtime
  // (note: does not include outfile)
  for(int i=0; i<GLOBAL_nfiles; i++){
    remove(GLOBAL_filenames[i]);
    free(GLOBAL_filenames[i]);
  }
  if(GLOBAL_filenames){
    free(GLOBAL_filenames);
    GLOBAL_filenames = NULL;
  }

  if(GLOBAL_all_leaves){
    free(GLOBAL_all_leaves);
    GLOBAL_all_leaves = NULL;
  }
  free_trie(GLOBAL_trie);
  GLOBAL_trie = NULL;

  if(GLOBAL_mergebuffers){
    for(int i=0; i<GLOBAL_nbuffers; i++){
      if(GLOBAL_mergebuffers[i]){
        free(GLOBAL_mergebuffers[i]);
        GLOBAL_mergebuffers[i] = NULL;
      }
    }
    if(GLOBAL_mergebuffers){
      free(GLOBAL_mergebuffers);
      GLOBAL_mergebuffers = NULL;
    }
  }

  if(GLOBAL_mergetree){
    LT_free(GLOBAL_mergetree);
    GLOBAL_mergetree = NULL;
  }

  if(GLOBAL_readedges){
    free(GLOBAL_readedges);
    GLOBAL_readedges = NULL;
  }
  // zero out all remaining constants
  GLOBAL_num_files = 0;
  GLOBAL_nfiles = 0;
  GLOBAL_cachectr = 0;
  return;
}

static inline char get_buffchar(char *buf, size_t bufsize, size_t *cur_ind,
                                size_t *remaining, FILE *stream){
  /*
   * buffered reading
   * if elements are in the buffer, returns the next buffered element
   * otherwise, refills the buffer with the next data from the file
   */

  if(*cur_ind == *remaining){
    *cur_ind = 0;
    *remaining = fread(buf, 1, bufsize, stream);
  }
  return *remaining ? buf[(*cur_ind)++] : 0;
}

static void truncate_file(const char* fname, size_t size){
  // not guaranteed to work, but not a big deal if it doesn't
  int retval = 0;
#ifdef HAVE_UNISTD_H
  retval = truncate(fname, size);
#else
  #ifdef WIN32
    FILE* f = fopen(fname, "rb+");
    int filehandler = _fileno(f);
    retval = _chsize_s(filehandler, size);
    fclose(f);
  #else
    return;
  #endif
#endif
  if(retval != 0) error("Failed to truncate file!");
  return;
}

/************************/
/* Comparison Functions */
/************************/

static int leaf_index_compar(const void *a, const void *b){
  l_uint aa = GLOBAL_leaf[*(l_uint*)a]->count;
  l_uint bb = GLOBAL_leaf[*(l_uint*)b]->count;
  return (aa > bb) - (aa < bb);
}

/*
 * Old comparator function in case this functionality is someday needed
static int edge_compar(const void *a, const void *b){
  // sort edges by INCREASING index, breaking ties by DECREASING weight
  const edge aa = *(const edge *)a;
  const edge bb = *(const edge *)b;
  int cmp = (aa.v > bb.v) - (aa.v < bb.v);
  // break ties by decreasing weight
  if(cmp == 0) {
    w_float w1, w2;
    l_uint v1, v2;
    decompressEdgeValue(aa.w, &v1, &w1);
    decompressEdgeValue(bb.w, &v2, &w2);
    w1 -= w2;
    cmp = (w1 < 0) - (w1 > 0);
    // break ties by increasing vertex 2 index
    if(cmp == 0) cmp = (v1 > v2) - (v1 < v2);
  }

  return cmp;
}
*/

static int fast_edge_compar(const void *a, const void *b){
  // skipping the decompression step to just directly sort
  // this sorts by increasing index 1 -> increasing index 2 -> undefined weight order
  const edge aa = *(const edge *)a;
  const edge bb = *(const edge *)b;
  int cmp = (aa.v > bb.v) - (aa.v < bb.v);
  if(cmp == 0) cmp = (aa.w > bb.w) - (aa.w < bb.w);
  return cmp;
}

/**************************/
/* File Mergesort Helpers */
/**************************/

static void kway_mergesort_file_inplace(const char* f1,
                          l_uint nlines,
                          size_t element_size,
                          l_uint block_size, int buf_size,
                          int num_bins, int output_size,
                          int (*compar)(const void *, const void *),
                          const int verbose){
  /*
   * In-place external merge sort.
   * Note that this scales quadratically in the worst case because we have
   *  to reorganize blocks as we go. Could do something smarter, but the
   *  literature is really complicated.
   *
   * Input variables:
   *  -         f1, f2: file names (f1: source, f2: temp)
   *  -         nlines: number of lines to sort in file
   *  -   element_size: size (in bytes) of a single element
   *  -     block_size: number of elements in each pre-sorted block
   *  -       buf_size: number of elements to store in each buffer
   *  -       num_bins: "k" in k-way merge. Set to power of 2 for best performance.
   *  -    output_size: number of elements to store in output buffer.
   *  -         compar: comparison function for comparing elements
   *  -        verbose: print status?
   */

  // file should already be sorted into x blocks of size block_size*element_size
  FILE *fileptr;

  LoserTree *mergetree = LT_alloc(num_bins, output_size, element_size, compar);
  GLOBAL_mergetree = mergetree;
  size_t cur_start, to_read;
  int empty_bin, LT_total_bins;
  l_uint nblocks, num_iter;
  double prev_progress, cur_progress;
  l_uint cur_fsize=0, max_fsize=0;

  if(num_bins > (nlines / block_size)){
    num_bins = nlines/block_size + !!(nlines % block_size);
  }
  if(verbose >= VERBOSE_BASIC) Rprintf("\tSorting with %d-way merge...\n", num_bins);

  l_uint nmax_iterations = 1, tmpniter=block_size;
  while(tmpniter*num_bins < nlines){
    tmpniter *= num_bins;
    nmax_iterations++;
  }
  tmpniter = 0;
  LT_total_bins = mergetree->nbins;

  // allocate space for the buffers
  void **buffers = safe_malloc(sizeof(void*)*num_bins);
  for(int i=0; i<num_bins; i++) buffers[i] = safe_malloc(buf_size*element_size);
  GLOBAL_mergebuffers = buffers;
  GLOBAL_nbuffers = num_bins;

  long int *offsets, *remaining;
  offsets = safe_calloc(LT_total_bins, sizeof(long int));
  remaining = safe_calloc(LT_total_bins, sizeof(long int));

  fileptr = safe_fopen(f1, "rb+");

  while(block_size < nlines){
    cur_progress = 0.0;
    prev_progress = 0.0;

    tmpniter++;
    if(verbose == VERBOSE_ALL){
      Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                              tmpniter, nmax_iterations, cur_progress);
      report_filesize(max_fsize);
      Rprintf(")  \r");
    } else if(verbose == VERBOSE_BASIC){
      Rprintf("\tIteration %" lu_fprint " of %" lu_fprint "\n", tmpniter, nmax_iterations);
    }
    R_CheckUserInterrupt();

    // number of blocks
    nblocks = nlines / block_size + !!(nlines % block_size);
    // number of iterations is nblocks / num bins
    num_iter = nblocks / num_bins + !!(nblocks % num_bins);

    cur_start = 0;
    for(l_uint iter=0; iter<num_iter; iter++){
      // run one k-way merge operation

      // first initialize offsets and number of remaining values
      for(int i=0; i<num_bins; i++){
        offsets[i] = cur_start;
        if(cur_start != nlines)
          cur_start += block_size;
        cur_start = cur_start > (nlines) ? nlines : cur_start;
        remaining[i] = cur_start - offsets[i];
      }

      // load data into the buffers and assign into tree
      for(int i=0; i<num_bins; i++){
        if(remaining[i]){
          fseek(fileptr, offsets[i]*element_size, SEEK_SET);
          to_read = remaining[i] > buf_size ? buf_size : remaining[i];
          safe_fread(buffers[i], element_size, to_read, fileptr);
          LT_fillBin(mergetree, i, to_read, buffers[i]);
          remaining[i] -= to_read;
          offsets[i] += to_read;
        }
      }

      // merge all the data
      LT_initGame(mergetree);
      while(mergetree->full_bins){
        // note cur_start is the first line of the *next* block
        empty_bin = LT_runInplaceFileGame(mergetree, cur_start, fileptr, remaining, &offsets);

        // refill the bin
        to_read = 0;
        if(remaining[empty_bin]){
          fseek(fileptr, offsets[empty_bin]*element_size, SEEK_SET);
          to_read = remaining[empty_bin] > buf_size ? buf_size : remaining[empty_bin];
          safe_fread(buffers[empty_bin], element_size, to_read, fileptr);
          remaining[empty_bin] -= to_read;
          offsets[empty_bin] += to_read;
        }

        cur_progress = ((double)mergetree->nwritten)/nlines*100.0;
        if(cur_progress - prev_progress > 0.5){
          prev_progress = cur_progress;
          if(verbose == VERBOSE_ALL){
            if(tmpniter == 1){
              cur_fsize = ftell(fileptr);
              if(cur_fsize > max_fsize) max_fsize = cur_fsize;
            }
            Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                              tmpniter, nmax_iterations, cur_progress);
            report_filesize(max_fsize);
            Rprintf(")  \r");
          }
          R_CheckUserInterrupt();
        }
        // note that the bin must be refilled even with no elements
        // so the tree moves on to the next step prior to next pop
        LT_refillBin(mergetree, empty_bin, to_read, buffers[empty_bin]);
      }

      // finally, call fdumpOutput to dump any remaining values
      // this doesn't need to call fdumpOutputInplace because there are no
      // values left to be saved from this block
      fseek(fileptr, (mergetree->nwritten)*element_size, SEEK_SET);
      LT_fdumpOutput(mergetree, fileptr);
      cur_progress = ((double)mergetree->nwritten)/nlines*100.0;

      if(cur_progress - prev_progress > 0.5 || cur_progress == 100){
          prev_progress = cur_progress;
          if(verbose == VERBOSE_ALL){
            if(tmpniter == 1){
              cur_fsize = ftell(fileptr);
              if(cur_fsize > max_fsize) max_fsize = cur_fsize;
            }
            Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                              tmpniter, nmax_iterations, cur_progress);
            report_filesize(max_fsize);
            Rprintf(")  \r");
          } else if(verbose == VERBOSE_BASIC && cur_progress == 100){
            Rprintf("\tIteration complete (used ");
            report_filesize(ftell(fileptr));
            Rprintf(")\n");
          }
        R_CheckUserInterrupt();
      }
    }

    mergetree->nwritten = 0;
    block_size *= num_bins;
    if(block_size >= nlines){
      cur_fsize = ftell(fileptr);
      if(cur_fsize > max_fsize) max_fsize = cur_fsize;
    }
    rewind(fileptr);
  }
  // cur_source will always be the file we just WROTE to here

  if(verbose >= VERBOSE_BASIC){
    Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                            tmpniter, nmax_iterations, 100.0);
    report_filesize(max_fsize);
    Rprintf(")  \n");
  }
  for(int i=0; i<num_bins; i++) free(buffers[i]);
  free(buffers);
  fclose_tracked(1);
  GLOBAL_mergebuffers = NULL;
  GLOBAL_nbuffers = 0;
  LT_free(mergetree);
  GLOBAL_mergetree = NULL;

  return;
}

static void kway_mergesort_file(const char* f1, const char* f2,
                          l_uint nlines,
                          size_t element_size,
                          l_uint block_size, int buf_size,
                          int num_bins, int output_size,
                          int (*compar)(const void *, const void *),
                          const int verbose){
  /*
   * Not-in-place external merge sort
   * Much faster than the in-place version, but takes double the disk space.
   *
   * Input variables:
   *  -       f1, f2: file names (f1: source, f2: temp)
   *  -       nlines: number of lines to sort in file
   *  - element_size: size (in bytes) of a single element
   *  -   block_size: number of elements in each pre-sorted block
   *  -     buf_size: number of elements to store in each buffer
   *  -     num_bins: "k" in k-way merge. Set to power of 2 for best performance.
   *  -  output_size: number of elements to store in output buffer.
   *  -       compar: comparison function for comparing elements
   *  -      verbose: print status?
   *
   */

  // file should already be sorted into x blocks of size block_size*element_size

  FILE *file1, *file2;

  LoserTree *mergetree = LT_alloc(num_bins, output_size, element_size, compar);
  GLOBAL_mergetree = mergetree;
  size_t cur_start, to_read;
  int empty_bin, LT_total_bins;
  l_uint nblocks, num_iter;
  double prev_progress, cur_progress;
  l_uint cur_fsize=0, max_fsize=0;

  if(num_bins > (nlines / block_size)){
    num_bins = nlines/block_size + !!(nlines % block_size);
  }
  if(verbose >= VERBOSE_BASIC) Rprintf("\tSorting with %d-way merge...\n", num_bins);

  l_uint nmax_iterations = 1, tmpniter=block_size;
  while(tmpniter*num_bins < nlines){
    tmpniter *= num_bins;
    nmax_iterations++;
  }
  tmpniter = 0;
  LT_total_bins = mergetree->nbins;

  // allocate space for the buffers
  void **buffers = safe_malloc(sizeof(void*)*num_bins);
  for(int i=0; i<num_bins; i++) buffers[i] = safe_malloc(buf_size*element_size);
  GLOBAL_mergebuffers = buffers;
  GLOBAL_nbuffers = num_bins;

  long int *offsets, *remaining;
  offsets = safe_calloc(LT_total_bins, sizeof(long int));
  remaining = safe_calloc(LT_total_bins, sizeof(long int));

  const char *cur_source = f1;
  const char *cur_target = f2;
  const char *tmp_swap_char;
  while(block_size < nlines){
    cur_progress = 0.0;
    prev_progress = 0.0;

    tmpniter++;
    if(verbose == VERBOSE_ALL){
      Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                        tmpniter, nmax_iterations, cur_progress);
      report_filesize(max_fsize);
      Rprintf(")  \r");
    }
    R_CheckUserInterrupt();

    file1 = safe_fopen(cur_source, "rb");
    file2 = safe_fopen(cur_target, "wb");

    // number of blocks
    nblocks = nlines / block_size + !!(nlines % block_size);
    // number of iterations is nblocks / num bins
    num_iter = nblocks / num_bins + !!(nblocks % num_bins);

    cur_start = 0;
    for(l_uint iter=0; iter<num_iter; iter++){
      // run one k-way merge operation

      // first initialize offsets and number of remaining values
      for(int i=0; i<num_bins; i++){
        offsets[i] = cur_start;
        if(cur_start != nlines-1)
          cur_start += block_size;
        cur_start = cur_start > (nlines) ? nlines : cur_start;
        remaining[i] = cur_start - offsets[i];
      }

      // load data into the buffers and assign into tree
      for(int i=0; i<num_bins; i++){
        if(remaining[i]){
          fseek(file1, offsets[i]*element_size, SEEK_SET);
          to_read = remaining[i] > buf_size ? buf_size : remaining[i];
          safe_fread(buffers[i], to_read, element_size, file1);
          LT_fillBin(mergetree, i, to_read, buffers[i]);
          remaining[i] -= to_read;
          offsets[i] += to_read;
        }
      }

      // merge all the data
      LT_initGame(mergetree);
      while(mergetree->full_bins){
        empty_bin = LT_runFileGame(mergetree, file2);
        // refill the bin
        to_read = 0;
        if(remaining[empty_bin]){
          fseek(file1, offsets[empty_bin]*element_size, SEEK_SET);
          to_read = remaining[empty_bin] > buf_size ? buf_size : remaining[empty_bin];
          safe_fread(buffers[empty_bin], to_read, element_size, file1);
          remaining[empty_bin] -= to_read;
          offsets[empty_bin] += to_read;
        }
        cur_progress = ((double)mergetree->nwritten)/nlines*100.0;
        if(cur_progress - prev_progress > 0.5){
          prev_progress = cur_progress;
          if(verbose == VERBOSE_ALL){
            if(tmpniter == 1){
              cur_fsize = ftell(file1) + ftell(file2);
              if(cur_fsize > max_fsize) max_fsize = cur_fsize;
            }
            Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                              tmpniter, nmax_iterations, cur_progress);
            report_filesize(max_fsize);
            Rprintf(")  \r");
          }
          R_CheckUserInterrupt();
        }
        // note that the bin must be refilled even with no elements
        // so the tree moves on to the next step prior to next pop
        LT_refillBin(mergetree, empty_bin, to_read, buffers[empty_bin]);
      }
      cur_progress = ((double)mergetree->nwritten)/nlines*100.0;
      if(cur_progress - prev_progress > 0.5 || cur_progress == 100){
        prev_progress = cur_progress;
        if(verbose == VERBOSE_ALL){
          if(tmpniter == 1){
            cur_fsize = ftell(file1) + ftell(file2);
            if(cur_fsize > max_fsize) max_fsize = cur_fsize;
          }
          Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                            tmpniter, nmax_iterations, cur_progress);
          report_filesize(max_fsize);
          Rprintf(")  \r");
        } else if(verbose == VERBOSE_BASIC && cur_progress == 100){
          Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                            tmpniter, nmax_iterations, cur_progress);
          report_filesize(max_fsize);
          Rprintf(")  \n");
        }
        R_CheckUserInterrupt();
      }

      // finally, call fdumpOutput to dump any remaining values
      LT_fdumpOutput(mergetree, file2);
    }
    mergetree->nwritten = 0;
    block_size *= num_bins;
    if(block_size >= nlines){
      cur_fsize = ftell(file1) + ftell(file2);
      if(cur_fsize > max_fsize) max_fsize = cur_fsize;
    }
    fclose_tracked(2);
    tmp_swap_char = cur_source;
    cur_source = cur_target;
    cur_target = tmp_swap_char;
  }
  // cur_source will always be the file we just WROTE to here

  if(verbose >= VERBOSE_BASIC){
    Rprintf("\tIteration %" lu_fprint " of %" lu_fprint " (%5.01f%% Complete, used ",
                      tmpniter, nmax_iterations, 100.0);
    report_filesize(max_fsize);
    Rprintf(")  \n");
  }
  for(int i=0; i<num_bins; i++) free(buffers[i]);
  free(buffers);

  GLOBAL_mergebuffers = NULL;
  GLOBAL_nbuffers = 0;
  LT_free(mergetree);
  GLOBAL_mergetree = NULL;

  remove(cur_target);
  if(f1 != cur_source){
    // if we used an odd number of iterations,
    // the final file is the temp file
    // so we have to "swap" them
    rename(cur_source, f1);
  }

  return;
}

static void split_sorted_file(const char* nfilename, const char* wfilename,
                        l_uint nlines, int v){
  /*
   * File to split the kway-sorted file into two other files
   *
   * After sorting, file is a bunch of {long, long} values
   * where the second long is a compressed value containing {index, weight}.
   * We need to discard the first long value of each pair and replace the file
   * with just the second long value from each pair. After that, we split the
   * weights into a wfilename and the index into nfilename.
   */

  double *valuebuffer = safe_malloc(L_SIZE * FILE_READ_CACHE_SIZE);
  edge *edgebuffer = safe_malloc(sizeof(edge)*FILE_READ_CACHE_SIZE);

  FILE *f_r, *f_w;
  f_r = safe_fopen(nfilename, "rb");
  f_w = safe_fopen(nfilename, "rb+"); // wb will delete the file

  l_uint cur_line = 0;
  l_uint nread = 0;

  // discard the first long in each edge struct, overwriting the file as we go
  while(cur_line != nlines){
    // read lines into buffer
    nread = nlines - cur_line;
    if(nread > FILE_READ_CACHE_SIZE) nread = FILE_READ_CACHE_SIZE;
    safe_fread(edgebuffer, sizeof(edge), nread, f_r);

    // copy second value of each edge into the double buffer
    for(int i=0; i<nread; i++) valuebuffer[i] = edgebuffer[i].w;

    // overwrite entries at front of file
    safe_fwrite(valuebuffer, L_SIZE, nread, f_w);
    cur_line += nread;
  }
  fclose_tracked(2);

  // reduce file size to reclaim the wasted space
  truncate_file(nfilename, L_SIZE*nlines);

  FILE *wf_w;
  f_r = safe_fopen(nfilename, "rb");
  f_w = safe_fopen(nfilename, "rb+");
  wf_w = safe_fopen(wfilename, "wb");

  free(edgebuffer);
  w_float *weights = safe_malloc(W_SIZE * FILE_READ_CACHE_SIZE);
  l_uint *indices = safe_malloc(L_SIZE * FILE_READ_CACHE_SIZE);
  cur_line = 0;
  nread = 0;

  double cur_progress, prev_progress = 0;

  // split the compressed {index,weight} across the two files
  while(cur_line != nlines){
    // read lines into current buffer
    nread = nlines - cur_line;
    if(nread > FILE_READ_CACHE_SIZE) nread = FILE_READ_CACHE_SIZE;
    safe_fread(valuebuffer, L_SIZE, nread, f_r);

    // update values
    for(int i=0; i<nread; i++)
      decompressEdgeValue(valuebuffer[i], &(indices[i]), &(weights[i]));

    safe_fwrite(indices, L_SIZE, nread, f_w);
    safe_fwrite(weights, W_SIZE, nread, wf_w);
    cur_line += nread;
    cur_progress = 100*((double)cur_line) / nlines;
    if(cur_progress - prev_progress > 0.5 || cur_progress == 100){
      if(v == VERBOSE_ALL) Rprintf("\tTidying up edges (%5.01f%% Complete)\r", cur_progress);
      prev_progress = cur_progress;
      R_CheckUserInterrupt();
    }
  }
  if(v == VERBOSE_ALL) Rprintf("\n");
  fclose_tracked(3);

  free(weights);
  free(indices);
  free(valuebuffer);

  return;
}

static void copy_weightsfile_sig(const char* dest, const char* src,
                                  l_uint num_edges, const double w){
  w_float *restrict wbuf = safe_malloc(FILE_READ_CACHE_SIZE*W_SIZE);
  FILE *fd = safe_fopen(dest, "wb");
  FILE *fs = safe_fopen(src, "rb");

  l_uint remaining=num_edges, to_read=0;

  while(remaining){
    to_read = remaining > FILE_READ_CACHE_SIZE ? FILE_READ_CACHE_SIZE : remaining;
    safe_fread(wbuf, W_SIZE, to_read, fs);
    for(l_uint i=0; i<to_read; i++)
      wbuf[i] = w < 0 ? 0 : sigmoid_transform(wbuf[i], w);
    safe_fwrite(wbuf, W_SIZE, to_read, fd);
    remaining -= to_read;
  }

  free(wbuf);
  fclose_tracked(2);
  return;
}

/******************/
/* Main Functions */
/******************/

static void reindex_trie_and_write_counts(prefix *trie, l_uint max_seen,
                                          int verbose, l_uint* nseen){
  /*
   * Vertices are 0-indexed
   * reindexing the leaves within a global array to find them without querying
   * the trie every time.
   */
  if(!trie) return;
  uint8_t bits_remaining = trie->count1 + trie->count2;
  uint8_t ctr = 0;
  l_uint cur_index;
  if(trie->bmap1 & 1){
    leaf *l = (leaf*)(trie->child_nodes[ctr++]);
    cur_index = l->index;
    GLOBAL_all_leaves[cur_index] = l;

    (*nseen)++;
    if(((*nseen)+1) % PRINT_COUNTER_MOD == 0){
      R_CheckUserInterrupt();
      if(verbose == VERBOSE_ALL) Rprintf("\tProcessed %" lu_fprint " vertices\r", *nseen);
    }
  }

  // hopefully this can be tail-call optimized by the compiler, but who knows
  while(ctr < bits_remaining)
    reindex_trie_and_write_counts(trie->child_nodes[ctr++], max_seen, verbose, nseen);
}

static void reset_trie_clusters(l_uint num_v){
  // same as reindex_trie_and_write_counts, but with no side effect writes
  // should already be initialized, just reset clusters to index+1 for all leaves
  leaf *l;
  for(l_uint i=0; i<num_v; i++){
    l = GLOBAL_all_leaves[i];
    l->count = l->index+1;
    l->dist = 0;
  }

  return;
}

static l_uint csr_compress_edgelist_trie(const char* edgefile, prefix *trie,
                                  FILE* neighbortable,
                                  const char sep, const char linesep,
                                  l_uint* num_v, int v,
                                  const int is_undirected,
                                  const int ignore_weights,
                                  int is_final_file){
  /*
   * This combines reading nodes and edges
   *
   * Note that the buffer of read edges is a global variable
   * (GLOBAL_readedges, GLOBAL_cachectr)
   * This is because we may have files that don't line up with the cache size,
   * leading to errors in the k-way mergesort later
   * (which expects sorted blocks of fixed size).
   * E.g. if we have two files of length 15 and cache size of 10, we'd write
   * four sorted blocks of size 10, 5, 10, 5, but k-way merge would expect 3
   * sorted of size 10, 10, 10. Using a global cache solves this issue by
   * allowing reads across files.
   */

  char *read_cache;
  char *restrict vname[2];
  char *restrict weight_buf;
  size_t rc_size = FILE_READ_CACHE_SIZE*8;

  l_uint nedges = 0;
  l_uint inds[2];
  float weight;

  int stringctr=0;
  l_uint print_counter = 0;
  size_t rcache_i = 0;
  size_t remaining = 0;

  // allocate space for all the stuff
  read_cache = safe_malloc(rc_size);
  vname[0] = safe_malloc(MAX_NODE_NAME_SIZE);
  vname[1] = safe_malloc(MAX_NODE_NAME_SIZE);
  if(!ignore_weights) weight_buf = safe_malloc(MAX_NODE_NAME_SIZE);

  FILE *edgelist = safe_fopen(edgefile, "rb");

  // note: weights/neighbors don't need to be initialized
  // see https://stackoverflow.com/questions/31642389/fseek-a-newly-created-file

  if(v >= VERBOSE_BASIC) Rprintf("\tReading file %s...\n", edgefile);

  char c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
  while(!feof(edgelist) || remaining){ // file isn't done OR buffer has data
    // read in the two vertex names
    for(int i=0; i<2; i++){
      stringctr = 0;
      while(c != sep && c != linesep){
        if(stringctr == MAX_NODE_NAME_SIZE)
          error("Incomplete entry read (suspect line %" lu_fprint ")", print_counter+1);
        vname[i][stringctr++] = c;
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
      }
      // short circuit to skip a length-0 name
      if(stringctr == 0){
        i--;
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
        continue;
      }
      vname[i][stringctr] = 0;

      // advance one past the separator if it isn't linesep
      // it would equal linesep if we don't have weights included
      if(c == sep)
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
    }
    stringctr = 0;
    if(!ignore_weights){
      // read in the weight
      memset(weight_buf, 0, MAX_NODE_NAME_SIZE);

      // need the double check in case user forgets a trailing \n
      while(c != linesep && (remaining || !feof(edgelist))){
        if(stringctr == MAX_NODE_NAME_SIZE)
          error("Incomplete entry read (suspect line %" lu_fprint ")", print_counter+1);
        weight_buf[stringctr++] = c;
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
      }
      weight = atof(weight_buf);
    } else {
      weight = 1.0;
      while(c != linesep && (remaining || !feof(edgelist)))
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
    }

    if(weight <= 0.0){
      // if weight is negative or zero, store the vertex but don't add an edge
      for(int i=0; i<2; i++)
        find_index_for_prefix_and_increment(vname[i], trie, num_v, 0);
    } else {
      // otherwise, add the edges to the cache
      for(int i=0; i<2; i++)
        inds[i] = find_index_for_prefix_and_increment(vname[i], trie, num_v, i + is_undirected);

      // add to cache, only adding the incoming edge if directed
      GLOBAL_readedges[GLOBAL_cachectr++] = compressEdgeValues(inds[1], inds[0], weight);
      if(is_undirected)
        GLOBAL_readedges[GLOBAL_cachectr++] = compressEdgeValues(inds[0], inds[1], weight);
    }

    // advance one past the separator
    c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
    if((GLOBAL_cachectr+is_undirected) >= FILE_READ_CACHE_SIZE || (is_final_file && feof(edgelist) && !remaining)){
      // sort the block and write it to the neighbors file
      qsort(GLOBAL_readedges, GLOBAL_cachectr, sizeof(edge), fast_edge_compar);
      nedges += safe_fwrite(GLOBAL_readedges, sizeof(edge), GLOBAL_cachectr, neighbortable);
      GLOBAL_cachectr = 0;
    }

    print_counter++;
    if(!(print_counter % PRINT_COUNTER_MOD)){
      if(v == VERBOSE_ALL) Rprintf("\t%" lu_fprint " lines read\r", print_counter);
      R_CheckUserInterrupt();
    }
  }
  if(v >= VERBOSE_BASIC) Rprintf("\t%" lu_fprint " lines read\n", print_counter);

  fclose_tracked(1);
  free(vname[0]);
  free(vname[1]);
  if(!ignore_weights) free(weight_buf);
  free(read_cache);

  return nedges;
}

static void add_remaining_to_queue(l_uint new_clust, leaf **neighbors,
                                    float *weights, l_uint nedge,
                                    ArrayQueue *queue){
  // called from update_node_cluster, adds all remaining leaves with different clusters to file
  l_uint ctr=0, tmp_ind, tmp_cl;
  int found;

  for(l_uint i=0; i<nedge; i++){
    tmp_cl = neighbors[i]->count;
    if(tmp_cl == new_clust) continue;
    tmp_ind = neighbors[i]->index;
    found = 0;
    for(l_uint j=0; j<ctr; j++){
      if(neighbors[j]->index == tmp_ind){
        found = 1;
        break;
      }
    }
    if(!found){
      neighbors[ctr++] = neighbors[i];
      if(ctr == MAX_EDGES_EXACT) break;
    }
  }

  for(l_uint j=0; j<ctr; j++)
    array_queue_insert(queue, neighbors[j]->index);

  return;
}

static void update_node_cluster(l_uint ind,
                          float self_loop_weight, float dist_pow,
                          FILE *weightsfile, FILE *neighborfile,
                          ArrayQueue *queue, float atten_param){
  /*
   * Determine number of edges using the table file (next - cur)
   * Inputs:
   *  -              ind: 0-indexed vertex id
   *  - self_loop_weight: self loop weight
   *  -      weightsfile: file pointer to weights for each edge
   *  -     neighborfile: file pointer to start nodes for each edge
   *  -            queue: queue of nodes to process
   *  -      atten_param: current attenuation parameter
   *
   */
  R_CheckUserInterrupt();
  l_uint start, end, num_edges;
  w_float *weights_arr;
  l_uint *indices;
  char *sufficient_weight;
  leaf **neighbors;
  leaf* original_node = GLOBAL_all_leaves[ind];
  l_uint original_cluster = original_node->count;

  // get the number of edges
  start = GLOBAL_all_leaves[ind]->edge_start;
  end = GLOBAL_all_leaves[ind+1]->edge_start;

  num_edges = end - start;
  // if it has no edges we can't do anything
  if(!num_edges) return;


  weights_arr = safe_malloc(W_SIZE*(num_edges));
  indices = safe_malloc(L_SIZE*(num_edges));
  neighbors = safe_malloc(sizeof(leaf*)*(num_edges));
  sufficient_weight = safe_malloc(num_edges);

  // read in the edges
  fseek(weightsfile, start*W_SIZE, SEEK_SET);
  fseek(neighborfile, start*L_SIZE, SEEK_SET);

  // these are the indexes and weights of the neighbors
  safe_fread(indices, L_SIZE, num_edges, neighborfile);
  safe_fread(weights_arr, W_SIZE, num_edges, weightsfile);

  // read in the clusters
  for(l_uint i=0; i<num_edges; i++){
    neighbors[i] = GLOBAL_all_leaves[indices[i]];
    sufficient_weight[i] = weights_arr[i] >= self_loop_weight;

    // attenuate edges, using adaptive scaling
    // see https://doi.org/10.1103/PhysRevE.83.036103, eqns 4-5
    // Note all weights are positive
    if(dist_pow != 1)
      weights_arr[i] *= 1-(atten_param * pow(neighbors[i]->dist, dist_pow));
    else
      weights_arr[i] *= 1-(atten_param * neighbors[i]->dist);
  }

  // sort both leaves and weights by assigned cluster
  // note that this just sorts the indexes vector, not the array itself
  for(l_uint i=0; i<num_edges; i++) indices[i] = i;
  GLOBAL_leaf = neighbors;
  qsort(indices, num_edges, L_SIZE, leaf_index_compar);


  // figure out the cluster to update to
  // TODO: shouldn't max_weight be DOUBLE_MIN ?
  double max_weight=MIN_DOUBLE, cur_weight=0, cur_error=0;
  l_uint max_clust=0, cur_clust=0;
  dist_uint new_dist = -1, min_dist = -1; // type defined in PrefixTrie.h
  char found_sufficient_weight = 0; // check if seen a weight at least self_loop_weight
  leaf* cur_neighbor;
  for(l_uint i=0; i<num_edges; i++){
    cur_neighbor = neighbors[indices[i]];
    if(cur_neighbor->count != cur_clust){
      // ensure we've found an edge at least larger than self_loop
      if(found_sufficient_weight && max_weight < cur_weight){
        max_weight = cur_weight;
        max_clust = cur_clust;
        new_dist = min_dist;
      }
      cur_clust = cur_neighbor->count;
      cur_weight = 0;
      min_dist = -1; // unsigned, so this becomes the maximum value
      cur_error = 0;
      found_sufficient_weight = 0;
    }

    kahan_accu(&cur_weight, &cur_error, weights_arr[indices[i]]);
    min_dist = min_dist > cur_neighbor->dist ? cur_neighbor->dist : min_dist;
    if(sufficient_weight[indices[i]])
      found_sufficient_weight = 1;
  }
  if(found_sufficient_weight && max_weight < cur_weight){
    max_weight = cur_weight;
    max_clust = cur_clust;
    new_dist = min_dist;
  }
  if(max_clust == 0){
    // max_clust == 0 case is for when the node has no edges,
    // which would skip the above loop and assign the node to 0 (invalid),
    // or if we never update it because all weights are negative
    // (due to attenuation) and none were originally larger than self-loop
    max_clust = original_cluster;
  }
  free(sufficient_weight);

  // have to actually write the new cluster and add changed nodes
  // only need to do this if it's changed, though
  if(max_clust != original_cluster){
    GLOBAL_verts_changed++;
    original_node->count = max_clust;
    // increment while handling overflow
    if(new_dist < DIST_UINT_MAX) new_dist++;
    original_node->dist = new_dist;
    add_remaining_to_queue(max_clust, neighbors, weights_arr, num_edges, queue);
  }

  free(weights_arr);
  free(neighbors);
  free(indices);
  return;
}

static void cluster_file(const char* weights_fname,
                          const char* neighbor_fname,
                          const l_uint num_v, const int max_iterations, const int v,
                          const float self_loop_weight,
                          const double atten_pow, const double dist_pow){
  // main runner function to cluster nodes
  FILE *weightsfile = safe_fopen(weights_fname, "rb+");
  FILE *neighborfile = safe_fopen(neighbor_fname, "rb");

  const char* progress[] = {"|o-----|", "|-o----|", "|--o---|", "|---o--|", "|----o-|", "|-----o|",
                            "|-----o|", "|----o-|", "|---o--|", "|--o---|", "|-o----|", "|o-----|"};
  const int progbarlen = 12;
  const int progbarcharlen = 8;
  const float MIN_UPDATE_PCT = 0.001;

  char progpreprint[progbarcharlen+1];
  memset(progpreprint, '\b', progbarcharlen);

  int print_counter=0, i=0;
  uint statusctr=0;
  l_uint iteration_length;
  float pct_complete = 0, prev_pct=0;

  // These are defined in PrefixTrie.h
  float atten_param = 0;

  // randomly initialize queue and ctr file

  if(v >= VERBOSE_BASIC) Rprintf("\tInitializing queues...");
  GLOBAL_queue = init_array_queue(num_v, max_iterations);
  iteration_length = GLOBAL_queue->length;
  if(v >= VERBOSE_BASIC) Rprintf("done.\n\tClustering network:\n");

  if(v == VERBOSE_ALL) Rprintf("\t0.0%% Complete %s", progress[++statusctr%progbarlen]);

  while(GLOBAL_queue->length){
    if(!iteration_length){
      atten_param = (float)GLOBAL_verts_changed / num_v;
      if(atten_pow == 0){
        atten_param = 0;
      } else if (atten_pow != 1){
        atten_param = pow(atten_param, atten_pow);
      }
      iteration_length = GLOBAL_queue->length;
      GLOBAL_verts_changed = 0;
    }
    print_counter++;
    if(!(print_counter % PROGRESS_COUNTER_MOD)){
      if(v == VERBOSE_ALL){
        pct_complete = ((float)(num_v - GLOBAL_queue->length)) / num_v;
        if(pct_complete != 1 && fabs(pct_complete - prev_pct) < MIN_UPDATE_PCT){
          pct_complete = prev_pct;
        } else {
          prev_pct = pct_complete;
        }
        Rprintf("\r\t%0.1f%% Complete %s", (pct_complete)*100, progress[++statusctr%progbarlen]);
      } else{
        R_CheckUserInterrupt();
      }
    }
    l_uint next_vert = array_queue_pop(GLOBAL_queue);
    iteration_length--;

    // will either be negative (because just removed from queue)
    // or zero (because seen max_iteration times)
    aq_int num_seen = -1*GLOBAL_queue->seen[next_vert];
    if(!num_seen) num_seen = max_iterations;
    update_node_cluster(next_vert, self_loop_weight, dist_pow,
                        weightsfile, neighborfile,
                        GLOBAL_queue, atten_param);
  }
  if(v >= VERBOSE_BASIC){
    if(v == VERBOSE_ALL) Rprintf("\r");
    if(max_iterations > 0)
      Rprintf("\t100%% Complete!                \n");
    else
      Rprintf("\tComplete! (%d total iterations)     \n", i+1);
  }
  free_array_queue(GLOBAL_queue);
  GLOBAL_queue = NULL;
  fclose_tracked(2);

  return;
}

static void resolve_cluster_consensus(const char* clusters,
                                      const char* weightsfile, const char* neighborfile,
                                      l_uint num_v, int num_runs){
  // overwrite the files with new info
  FILE *clusts = safe_fopen(clusters, "rb");

  l_uint *read_clusts = safe_malloc(L_SIZE*num_v);
  l_uint *counts = safe_malloc(L_SIZE*num_v);
  l_uint *tmp_space = safe_malloc(L_SIZE*num_v);
  l_uint tmp, total_edges;

  memset(counts, 0, L_SIZE*num_v);

  for(int i=0; i<num_runs; i++){
    // read in clusters
    safe_fread(read_clusts, L_SIZE, num_v, clusts);

    // tabulate clusters
    memset(tmp_space, 0, L_SIZE*num_v);
    for(l_uint j=0; j<num_v; j++)
      tmp_space[read_clusts[j]-1]++;

    for(l_uint j=0; j<num_v; j++)
      // number of elements in this cluster must be at least 1 if this node is in it
      counts[j] += tmp_space[read_clusts[j]-1] - 1;
  }

  // convert to cumulative counts
  tmp = 0;
  for(l_uint i=0; i<num_v; i++){
    GLOBAL_all_leaves[i]->edge_start = tmp;
    tmp += counts[i];
  }
  total_edges = tmp;

  // write cumulative counts to file
  l_uint cur_clust, tmp_ind, num_neighbors;
  FILE *neighbors = safe_fopen(neighborfile, "wb+");

  // having fwrite problems, I'm going to just set up the file in advance
  tmp_ind = total_edges;
  for(l_uint i=0; i<num_v; i++)
    tmp_space[i] = 0;
  while(tmp_ind){
    tmp = tmp_ind > num_v ? num_v : tmp_ind;
    safe_fwrite(tmp_space, L_SIZE, tmp, neighbors);
    tmp_ind -= tmp;
  }

  for(int i=0; i<num_runs; i++){
    // read in all clusters
    safe_fread(read_clusts, L_SIZE, num_v, clusts);

    // tabulate clusters
    for(l_uint j=0; j<num_v; j++){
      if(!read_clusts[j]) continue;

      tmp = 1;
      tmp_space[0] = j;
      cur_clust = read_clusts[j];
      read_clusts[j] = 0;
      for(l_uint k=j+1; k<num_v; k++){
        // can only be later nodes, if it were earlier we would've already caught it
        if(read_clusts[k] == cur_clust){
          tmp_space[tmp++] = k;
          read_clusts[k] = 0;
        }
      }

      num_neighbors = tmp-1;
      if(!num_neighbors) continue;
      // now all the elements in the same cluster are in tmp_space[0:tmp-1]
      for(l_uint k=0; k<tmp; k++){
        // swap the current element to write to the beginning
        cur_clust = tmp_space[0];
        tmp_space[0] = tmp_space[k];
        tmp_space[k] = cur_clust;

        // write all the neighbors to the file
        cur_clust = tmp_space[0];
        tmp_ind = GLOBAL_all_leaves[cur_clust]->edge_start;

        // decrement counts first
        counts[cur_clust] -= num_neighbors;

        // then use it as index
        tmp_ind += counts[cur_clust];

        // write the weights later
        fseek(neighbors, tmp_ind*L_SIZE, SEEK_SET);
        safe_fwrite(&(tmp_space[1]), L_SIZE, num_neighbors, neighbors);
      }
    }
  }
  fclose_tracked(1);

  free(tmp_space);
  free(read_clusts);
  free(counts);

  // write the weights, should think of a better way to do ignore_weights
  FILE *weights = safe_fopen(weightsfile, "wb");
  w_float *w = safe_malloc(W_SIZE*FILE_READ_CACHE_SIZE);
  for(int i=0; i<FILE_READ_CACHE_SIZE; i++)
    w[i] = 1.0;
  while(total_edges){
    tmp = FILE_READ_CACHE_SIZE > total_edges ? total_edges : FILE_READ_CACHE_SIZE;
    safe_fwrite(w, W_SIZE, tmp, weights);
    total_edges -= tmp;
  }
  free(w);
  fclose_tracked(2);

  return;
}

static void consensus_cluster_oom(const char* weightsfile,
                            const char* neighborfile, const char* dir,
                            const l_uint num_v, const int num_iter, const int v,
                            const float self_loop_weight,
                            const double atten_pow, const double dist_pow,
                            const double* consensus_weights, const int consensus_len){

  /*
   * Inputs:
   *  - weightsfile: weights of each edge
   *  - neighborfile: neighbors of each edge
   *
   * Need to do the following:
   *  - copy weightsfile
   *  - apply any weights transformations (can combine with previous step)
   *  - re-initialize clusters
   *  - run cluster_file using weights file (clusters stored in trie)
   *  - store clusters somewhere
   */
  const char* transformedweights = create_filename(dir, CONSENSUS_CSRCOPY1);
  const char* tmpclusterfile = create_filename(dir, CONSENSUS_CLUSTER);

  FILE *dummyclust;

  // need to get the total number of edges
  l_uint num_edges = GLOBAL_all_leaves[num_v]->edge_start;

  l_uint *clusters = safe_malloc(L_SIZE * num_v);
  leaf *tmpleaf;

  dummyclust = safe_fopen(tmpclusterfile, "wb");
  // now we run clustering over consensus_len times
  for(int i=0; i<consensus_len; i++){
    if(v >= VERBOSE_BASIC) Rprintf("Iteration %d of %d:\n", i+1, consensus_len);

    // modify weights according to sigmoid transformation
    if(v >= VERBOSE_BASIC) Rprintf("\tTransforming edge weights...\n");
    copy_weightsfile_sig(transformedweights, weightsfile, num_edges, consensus_weights[i]);

    // reset cluster values
    reset_trie_clusters(num_v);

    // cluster with transformed weights
    cluster_file(transformedweights, neighborfile,
                  num_v, num_iter, v, self_loop_weight, atten_pow, dist_pow);

    if(v >= VERBOSE_BASIC) Rprintf("\tRecording results...\n");
    for(l_uint i=0; i<num_v; i++){
      tmpleaf = GLOBAL_all_leaves[i];
      clusters[tmpleaf->index] = tmpleaf->count;
    }

    // this should be fine assuming 64-bit system
    // R isn't supported on 32-bit machines, so we know it'll be large enough
    safe_fwrite(clusters, L_SIZE, num_v, dummyclust);
  }
  fclose_tracked(1);
  free(clusters);

  // now we can just destroy the other files
  if(v >= VERBOSE_BASIC) Rprintf("Reconciling runs...\n");
  resolve_cluster_consensus(tmpclusterfile, weightsfile, neighborfile,
                            num_v, consensus_len);

  if(v >= VERBOSE_BASIC) Rprintf("Clustering on consensus data...\n");
  reset_trie_clusters(num_v);
  cluster_file(weightsfile, neighborfile, num_v, num_iter, v, self_loop_weight, atten_pow, dist_pow);

  return;
}

static l_uint write_output_clusters_trie(FILE *outfile, prefix *trie, l_uint *clust_mapping,
                                char *s, int cur_pos, char *write_buf, const size_t num_bytes,
                                const char *seps, l_uint num_v, int verbose){
  // traverse the tree to rebuild node names
  // get the clusters out of them and then write them
  // can optimize slightly by storing prefixes up to this point during a DFS
  if(!num_v) return 0;
  uint8_t bits_remaining = trie->count1 + trie->count2;
  uint8_t ctr = 0;
  uint8_t current_bit = 0;
  uint64_t bitmap = trie->bmap1;
  while(ctr < trie->count1){
    // iterate over first bitmap
    if(bitmap & 1){
      if(current_bit == 0){
        // read in the cluster (stored in leaf node at count)
        l_uint tmpval = ((leaf *)(trie->child_nodes[0]))->count;

        tmpval--; // decrement to make them 0-indexed
        if(!clust_mapping[tmpval]){
          clust_mapping[tmpval] = CLUST_MAP_CTR++;
        }
        tmpval = clust_mapping[tmpval];

        // prepare data for output
        s[cur_pos] = 0;
        snprintf(write_buf, num_bytes, "%s%c%" lu_fprint "%c", s, seps[0], tmpval, seps[1]);

        safe_fwrite(write_buf, 1, strlen(write_buf), outfile);
        num_v--;
        if(num_v % PROGRESS_COUNTER_MOD == 0){
          if(verbose == VERBOSE_ALL)
            Rprintf("\r\tNodes remaining: %" lu_fprint "                    ", num_v);
          else
            R_CheckUserInterrupt();
        }
      } else {
        // set character and recur
        s[cur_pos] = current_bit + CHAR_OFFSET;
        num_v = write_output_clusters_trie(outfile, trie->child_nodes[ctr], clust_mapping,
                                    s, cur_pos+1, write_buf, num_bytes, seps, num_v, verbose);
      }
      ctr++;
    }
    bitmap >>= 1;
    current_bit++;
  }

  bitmap = trie->bmap2;
  current_bit = 0;
  while(ctr < bits_remaining){
    // iterate over the second bitmap
    if(bitmap&1){
      // don't have to safecheck the bit==0 case
      s[cur_pos] = current_bit + MIN_VALUE_BMAP2;
      num_v = write_output_clusters_trie(outfile, trie->child_nodes[ctr], clust_mapping,
                                    s, cur_pos+1, write_buf, num_bytes, seps, num_v, verbose);
      ctr++;
    }
    bitmap >>= 1;
    current_bit++;
  }

  return num_v;
}

/***********/
/* TESTING */
/***********/

/*
 * These functions must be included here because most ExoLabel methods are static
 * These methods should allow validation of everything in here from R.
 */

/*
 * DISJOINT SETS
 */

static inline l_uint get_dsu_rep(l_uint* reps, l_uint v){
  if(reps[v] != v) reps[v] = get_dsu_rep(reps, reps[v]);
  return reps[v];
}

static void merge_dsu_sets(l_uint *reps, l_uint *sizes, int v1, int v2){
  l_uint key1 = get_dsu_rep(reps, v1);
  l_uint key2 = get_dsu_rep(reps, v2);
  if(key1 == key2) return;
  if(sizes[key1] > sizes[key2]){
    l_uint tmp = key1;
    key1 = key2;
    key2 = tmp;
  }
  reps[key1] = key2;
  sizes[key2] += sizes[key1];
  return;
}

static SEXP find_disjoint_sets(l_uint num_v, const double cutoff,
                              const char* edges_file, const char* weights_file){
  l_uint* sets = safe_malloc(sizeof(l_uint) * num_v);
  l_uint* sizes = safe_malloc(sizeof(l_uint) * num_v);
  FILE *efile = safe_fopen(edges_file, "rb");
  FILE *wfile = safe_fopen(weights_file, "rb");

  for(l_uint i = 0; i < num_v; i++){
    sets[i] = i;
    sizes[i] = 1;
  }

  l_uint start, nedges, buf_size = 0;
  l_uint *indices = NULL;
  w_float *weights = NULL;
  for(l_uint i=0; i<num_v; i++){
    start = GLOBAL_all_leaves[i]->edge_start;
    nedges = GLOBAL_all_leaves[i+1]->edge_start - start;
    if(nedges == 0) continue;
    if(buf_size < nedges){
      indices = safe_realloc(indices, nedges*L_SIZE);
      weights = safe_realloc(weights, nedges*W_SIZE);
      buf_size = nedges;
    }
    fseek(efile, start*L_SIZE, SEEK_SET);
    fseek(wfile, start*W_SIZE, SEEK_SET);
    safe_fread(indices, L_SIZE, nedges, efile);
    safe_fread(weights, W_SIZE, nedges, wfile);
    for(l_uint j=0; j<nedges; j++)
      if(weights[j] >= cutoff)
        merge_dsu_sets(sets, sizes, i, indices[j]);
  }

  if(indices) free(indices);
  if(weights) free(weights);

  // this will compress all paths to their final value
  l_uint num_disjoint_sets = 0;
  for(l_uint i=0; i<num_v; i++){
    if(get_dsu_rep(sets, i) != i)
      sizes[i] = 0;
    else
      num_disjoint_sets++;
  }
  SEXP groups = PROTECT(allocVector(REALSXP, num_disjoint_sets));
  double* to_return = REAL(groups);
  for(l_uint i=0; i<num_v; i++)
    if(sizes[i])
      to_return[--num_disjoint_sets] = (double)sizes[i];

  free(sets);
  free(sizes);
  fclose_tracked(2);
  return groups;
}

static SEXP validate_read_weights(const char* neighbor, l_uint expected_edges, int is_undirected){
  FILE *f = safe_fopen(neighbor, "rb");
  SEXP all_w = PROTECT(allocVector(REALSXP, expected_edges));
  double* result = REAL(all_w);
  edge* buf = safe_malloc(sizeof(edge)*(FILE_READ_CACHE_SIZE));
  l_uint v;
  w_float w;
  l_uint nread = FILE_READ_CACHE_SIZE;
  l_uint total_read = 0;
  while(nread == FILE_READ_CACHE_SIZE){
    nread = fread(buf, sizeof(edge), FILE_READ_CACHE_SIZE, f);
    for(int i=0; i<nread; i++){
      decompressEdgeValue(buf[i].w, &v, &w);
      result[total_read++] = w;
    }
  }
  fclose_tracked(1);
  free(buf);

  return all_w;
}

static SEXP validate_node_degree(l_uint num_v){
  SEXP all_degrees = PROTECT(allocVector(REALSXP, num_v));
  double* degs = REAL(all_degrees);
  l_uint tmp;
  for(l_uint i=0; i<num_v; i++){
    tmp = GLOBAL_all_leaves[i+1]->edge_start - GLOBAL_all_leaves[i]->edge_start;
    degs[i] = (double)tmp;
  }
  return all_degrees;
}

static void validate_edgefile_is_sorted(const char* neighbor, const l_uint expected_edges){
  // this function validates file contents after edges are sorted but before split
  FILE *f = safe_fopen(neighbor, "rb");
  edge* buf = safe_malloc(sizeof(edge)*(FILE_READ_CACHE_SIZE+1));
  buf[0].v = 0;

  l_uint nread = FILE_READ_CACHE_SIZE;
  l_uint total_read = 0;
  while(nread == FILE_READ_CACHE_SIZE){
    if(total_read > 0)
      buf[0].v = buf[FILE_READ_CACHE_SIZE].v;
    nread = fread(buf+1, sizeof(edge), FILE_READ_CACHE_SIZE, f);
    for(int i=1; i<=nread; i++){
      if(buf[i].v < buf[i-1].v) error("Edgefile is improperly sorted!\n");
    }
    total_read += nread;
  }
  fclose_tracked(1);
  free(buf);
  if(total_read != expected_edges){
    error("Wrong number of edges! Expected %" lu_fprint ", read %" lu_fprint ".\n",
      expected_edges, total_read);
  }
  return;
}


SEXP R_LPOOM_cluster(SEXP FILENAME, SEXP NUM_EFILES, // files
                    SEXP OUTDIR, SEXP OUTFILES,  // more files
                    SEXP SEPS, SEXP ITER, SEXP VERBOSE, // control flow
                    SEXP IS_UNDIRECTED, SEXP ADD_SELF_LOOPS, // optional adjustments
                    SEXP IGNORE_WEIGHTS, SEXP CONSENSUS_WEIGHTS,
                    SEXP SORT_INPLACE, SEXP ATTEN_POWER, SEXP DIST_POWER,
                    SEXP TESTING){
  /*
   * I always forget how to handle R strings so I'm going to record it here
   * R character vectors are STRSXPs, which is the same as a list (VECSXP)
   * Each entry in a STRSXP is a CHARSXP, accessed with STRING_ELT(str, ind)
   * You can get a `const char*` from a CHARSXP with CHAR()
   * You can also re-encode with `const char* Rf_translateCharUTF8()`
   */

  /*
   * Input explanation:
   *     FILENAME: file of edges in format `v1 v2 w`
   *       OUTDIR: directory to store temporary files
   */

  // initialize global variables
  GLOBAL_nfiles = 0;
  GLOBAL_filenames = safe_malloc(sizeof(char*) * 10);
  GLOBAL_ftracker = safe_malloc(sizeof(FILE*) * 20);
  GLOBAL_trie = initialize_trie();
  GLOBAL_all_leaves = NULL;
  GLOBAL_mergebuffers = NULL;
  GLOBAL_nbuffers = 0;
  GLOBAL_mergetree = NULL;
  GLOBAL_readedges = NULL;
  GLOBAL_cachectr = 0;

  // main files
  const char* dir = CHAR(STRING_ELT(OUTDIR, 0));
  const char* weightsfile = create_filename(dir, "weights.bin");
  const char* neighborfile = create_filename(dir, "neighbors.bin");
  const int num_ofiles = LENGTH(OUTFILES);
  const char* outfile;
  const char* edgefile;

  // required parameters
  const char* seps = CHAR(STRING_ELT(SEPS, 0));
  const int num_edgefiles = INTEGER(NUM_EFILES)[0];
  int* num_iter = INTEGER(ITER);
  aq_int base_iter = 0;
  const int verbose = INTEGER(VERBOSE)[0];
  l_uint num_v = 0;

  // optional parameters
  const int is_undirected = LOGICAL(IS_UNDIRECTED)[0];
  const double* self_loop_weights = REAL(ADD_SELF_LOOPS);
  const int ignore_weights = LOGICAL(IGNORE_WEIGHTS)[0];
  const int use_inplace_sort = LOGICAL(SORT_INPLACE)[0];
  const double* atten_power = REAL(ATTEN_POWER);
  const double* dist_power = REAL(DIST_POWER);
  const int debug = LOGICAL(TESTING)[0];

  // consensus stuff
  const int consensus_len = length(CONSENSUS_WEIGHTS);
  const double* consensus_w = REAL(CONSENSUS_WEIGHTS);

  // timing
  time_t time1, time2;

  time1 = time(NULL);
  if(verbose >= VERBOSE_BASIC) Rprintf("Reading in edges...\n");
  l_uint num_edges = 0;
  FILE *neighbortable = safe_fopen(neighborfile, "ab");
  GLOBAL_readedges = safe_malloc(FILE_READ_CACHE_SIZE * sizeof(edge));
  GLOBAL_cachectr = 0;
  for(int i=0; i<num_edgefiles; i++){
    edgefile = CHAR(STRING_ELT(FILENAME, i));
    num_edges += csr_compress_edgelist_trie(edgefile, GLOBAL_trie,
                                              neighbortable,
                                              seps[0], seps[1],
                                              &num_v, verbose,
                                              is_undirected,
                                              ignore_weights,
                                              i == num_edgefiles-1);
  }
  fclose_tracked(1);
  if(verbose >= VERBOSE_BASIC)
    print_graph_stats(num_v, num_edges);

  // allocate space for leaf counters (one extra space to hold final value of cumulative counts)
  GLOBAL_all_leaves = safe_malloc(sizeof(leaf *) * (num_v+1));

  // next, reformat the file to get final counts for each node
  if(verbose >= VERBOSE_BASIC) Rprintf("Tidying up internal tables...\n");
  l_uint print_val = 0;
  reindex_trie_and_write_counts(GLOBAL_trie, 0, verbose, &print_val);
  // final print in reindex_trie_and_write_counts is \r, need a newline
  if(verbose >= VERBOSE_ALL) Rprintf("\n");

  l_uint max_degree = 0;
  // change edge_start values to cumulative counts
  GLOBAL_all_leaves[num_v] = malloc(sizeof(leaf));
  GLOBAL_all_leaves[num_v]->count = 0;
  GLOBAL_all_leaves[num_v]->edge_start = 0;
  l_uint running_sum = 0;
  for(l_uint i=0; i<=num_v; i++){
    leaf* tmp_leaf = GLOBAL_all_leaves[i];
    if(tmp_leaf->count > max_degree) max_degree = tmp_leaf->count;
    if(running_sum > num_edges)
      error("Offsets incorrect: node %" lu_fprint
            " has offset %" lu_fprint "(only %"
            lu_fprint " edges total!)\n",
            i, running_sum, num_edges);
    tmp_leaf->edge_start = running_sum;
    running_sum += tmp_leaf->count;
    // trie clusters are reset later, no need to reset them here
  }
  if(verbose >= VERBOSE_ALL)
    Rprintf("\tMaximum node degree is %" lu_fprint "!\n", max_degree);

  SEXP debug_weights, debug_degrees;
  if(debug){
    if(verbose >= VERBOSE_BASIC) Rprintf("Recording node degrees...");
    debug_degrees = validate_node_degree(num_v);
    if(verbose >= VERBOSE_BASIC) Rprintf("done!\n");
    if(verbose >= VERBOSE_BASIC) Rprintf("Recording edge weights...");
    debug_weights = validate_read_weights(neighborfile, GLOBAL_all_leaves[num_v]->edge_start, is_undirected);
    if(verbose >= VERBOSE_BASIC) Rprintf("done!\n");
  }

  // get base number of iterations
  max_degree = (l_uint)(sqrt((double)max_degree));
  size_t num_bits = sizeof(aq_int) * 8 - 1; // signed, so we have one less bit to work with
  base_iter = ((aq_int)(1)) << num_bits;
  base_iter = base_iter > max_degree ? max_degree : base_iter;
  base_iter = base_iter < 5 ? 5 : base_iter; // minimum of 5 iterations per node

  // sort the file and split it into neighbor and weight
  if(FILE_READ_CACHE_SIZE < num_edges){
    if(use_inplace_sort){
      kway_mergesort_file_inplace(neighborfile,
                                  num_edges, sizeof(edge),
                                  FILE_READ_CACHE_SIZE,
                                  MERGE_INPUT_SIZE,
                                  MAX_BINS_FOR_MERGE, // bins to merge with
                                  MERGE_OUTPUT_SIZE, // size of output buffer
                                  fast_edge_compar, verbose);
    } else {
      kway_mergesort_file(neighborfile, weightsfile,
                                  num_edges, sizeof(edge),
                                  FILE_READ_CACHE_SIZE,
                                  MERGE_INPUT_SIZE,
                                  MAX_BINS_FOR_MERGE, // bins to merge with
                                  MERGE_OUTPUT_SIZE, // size of output buffer
                                  fast_edge_compar, verbose);
    }
  }

  if(debug){
    if(verbose >= VERBOSE_BASIC) Rprintf("Validating edges are sorted correctly...");
    validate_edgefile_is_sorted(neighborfile,  GLOBAL_all_leaves[num_v]->edge_start);
    if(verbose >= VERBOSE_BASIC) Rprintf("passed!\n");
  }

  split_sorted_file(neighborfile, weightsfile, num_edges, verbose);
  time2 = time(NULL);

  if(verbose >= VERBOSE_BASIC) report_time(time1, time2, "\t");

  SEXP debug_disjoint_sizes;
  if(debug){
    if(verbose >= VERBOSE_BASIC) Rprintf("Recording disjoint sets...");
    debug_disjoint_sizes = find_disjoint_sets(num_v, self_loop_weights[0],
                                            neighborfile, weightsfile);
    if(verbose >= VERBOSE_BASIC) Rprintf("done!\n");
  }

  char vert_name_holder[MAX_NODE_NAME_SIZE];
  char write_buffer[PATH_MAX];
  for(int i=0; i<num_ofiles; i++){
    if(verbose >= VERBOSE_BASIC) Rprintf("Clustering...\n");
    reset_trie_clusters(num_v);
    if(num_iter[i] == 0){
      if(verbose >= VERBOSE_BASIC)
        Rprintf("\tAutomatically setting iterations to %d (was 0)\n", base_iter);
      num_iter[i] = base_iter;
    }
    time1 = time(NULL);
    if(consensus_len){
      consensus_cluster_oom(weightsfile, neighborfile, dir, num_v,
                            num_iter[i], verbose, self_loop_weights[i],
                            atten_power[i], dist_power[i],
                            consensus_w, consensus_len);

    } else {
      cluster_file(weightsfile, neighborfile, num_v, num_iter[i], verbose,
                    self_loop_weights[i], atten_power[i], dist_power[i]);
    }
    time2 = time(NULL);
    if(verbose >= VERBOSE_BASIC) report_time(time1, time2, "\t");

    // have to allocate resources for writing out
    outfile = CHAR(STRING_ELT(OUTFILES, i));
    FILE *results = safe_fopen(outfile, "wb");
    if(!results) error("Failed to open output file.");
    l_uint *clust_mapping = safe_calloc(num_v, L_SIZE);
    CLUST_MAP_CTR = 1;
    if(verbose == VERBOSE_ALL) Rprintf("Writing clusters to file...\n\tVertices remaining: %" lu_fprint "", num_v);
    else if(verbose == VERBOSE_BASIC) Rprintf("Writing clusters to file...\n");
    write_output_clusters_trie(results, GLOBAL_trie, clust_mapping, vert_name_holder,
                                0, write_buffer, (size_t)PATH_MAX, seps, num_v, verbose);
    if(verbose == VERBOSE_ALL) Rprintf("\r\tNodes remaining: Done!               \n");
    else if(verbose == VERBOSE_BASIC) Rprintf("Execution completed!\n");
    free(clust_mapping);
    fclose_tracked(1);
  }

  remove(weightsfile);
  remove(neighborfile);
  cleanup_ondisklp_global_values();


  SEXP RETURN_VAL = PROTECT(allocVector(REALSXP, 2));
  REAL(RETURN_VAL)[0] = num_v;
  REAL(RETURN_VAL)[1] = num_edges;

  if(debug){
    SEXP LIST_VAL = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(LIST_VAL, 0, RETURN_VAL);
    SET_VECTOR_ELT(LIST_VAL, 1, debug_weights);
    SET_VECTOR_ELT(LIST_VAL, 2, debug_degrees);
    SET_VECTOR_ELT(LIST_VAL, 3, debug_disjoint_sizes);
    UNPROTECT(5);
    return LIST_VAL;
  }
  UNPROTECT(1);
  return RETURN_VAL;
}
