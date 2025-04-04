/*
 * Out of memory clustering with fast label propagation
 * Author: Aidan Lakshman
 *
 * This set of functions creates 4 files:
 *  -         csr: csr-compressed graph structure. Contains n+1 uint64_t values
 *                 corresponding to vertex indices, where the k'th value denotes
 *                 the start position in the other files for vertex k. For
 *                 example, if the first two values are 0 100, then the outgoing
 *                 edges from the first vertex are entries 0-99.
 *  -   neighbors: destination for each edge, indexed by `csr`.
 *  -     weights: weights for each edge, indexed by `csr`.
 *  - outfile.tsv: .tsv file returned to R, contains two tab-separated columns
 *                 (vertex name, cluster)
 *
 * Additional Notes:
 *  - sizeof(char) is guaranteed to be 1 (see https://en.wikipedia.org/wiki/Sizeof).
 *    1 is used instead to simplify code somewhat.
 *  - fwrite is thread safe but gives no performance benefit
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
 *  - switch back to uint16_t for node name length
 *
 * Optimization Ideas:
 *  - Can we speed up reading edges more?
 *    - add a cache of "recently used vertices" as [c1,c2,...,\0,ptr]
 *    - is it possible to order the nodes by what is most often referenced?
 *      ...seems like more trouble than its worth
 *    - reduce number of block interchanges in the in-place mergesort
 *  - Can we optimize RAM consumption?
 *    - might be worthwhile to compress the trie structure *after* its read in
 *  - Can we offload anything else to the trie?
 *      - could store offsets during clustering step so we don't have to query `csr`
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


// max size of a vertex name (char array will have 2 extra spaces for terminator and flag)
#define MAX_NODE_NAME_SIZE 254
// holds pointers char* of size MAX_NODE_NAME_SIZE, 4096 is 1MB
#define NODE_NAME_CACHE_SIZE 40960
// number of entries, so total consumption often multiplied by 8
#define FILE_READ_CACHE_SIZE 8192*4
#define CLUSTER_MIN_WEIGHT 0.01

/*************/
/* Constants */
/*************/

static const int L_SIZE = sizeof(l_uint);
static const int W_SIZE = sizeof(w_float);
static const char CONSENSUS_CSRCOPY1[] = "tmpcsr1";
static const char CONSENSUS_CLUSTER[] = "tmpclust";
static const int BITS_FOR_WEIGHT = 10;

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
 *  - Mergesort bin space is multiplied by sizeof(edge) (16 bytes)
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
  l_uint ctr1;
  l_uint ctr2;
} double_lu;

typedef struct {
  strlen_uint strlength;
  char s[MAX_NODE_NAME_SIZE];
  h_uint hash;
  l_uint count;
} msort_vertex_line;

typedef struct ll {
  l_uint id;
  float w;
  struct ll* next;
} ll;

typedef struct ll2 {
  float w;
  struct ll2* next;
} ll2;

typedef struct {
  strlen_uint len;
  h_uint hash;
  l_uint index;
} iline;

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
static l_uint GLOBAL_verts_remaining = 0;
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
static double GLOBAL_max_weight = 1.0;
static l_uint GLOBAL_verts_changed = 0;

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
  if(days) Rprintf("%d days, ", days);
  if(hours) Rprintf("%d hrs, ", hours);
  if(mins) Rprintf("%d mins, ", mins);
  Rprintf("%d secs\n", secs);
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

/************************/
/* Arithmetic Functions */
/************************/

static void kahan_accu(double *cur_sum, double *cur_err, double new_val){
  // this is an accumulator that will use a kahan sum to reduce floating point accuracy
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
  edge e;
  e.v = v1;
  l_uint v2_comp = v2;
  v2_comp <<= BITS_FOR_WEIGHT;

  // normalize to 0-1 scale
  weight = weight / GLOBAL_max_weight;

  // quantize to number of allowable bits
  l_uint max_weight_bits = (1 << BITS_FOR_WEIGHT) - 1;
  l_uint w_comp = (l_uint) floor(weight * max_weight_bits);

  e.w = v2_comp | w_comp;
  return e;
}

static void decompressEdgeValue(l_uint compressed, l_uint *v2, w_float *w){
  l_uint max_weight_bits = (1 << BITS_FOR_WEIGHT) - 1;
  l_uint w_comp = compressed & max_weight_bits;

  // two assignments because I'd like to minimize loss of precision
  double w_temp = ((double)w_comp) / max_weight_bits;
  w_temp *= GLOBAL_max_weight;

  //double w_temp = (((double)w_comp) * 2 * GLOBAL_max_weight) / max_weight_bits;
  *w = (w_float)w_temp;
  *v2 = compressed >> BITS_FOR_WEIGHT;

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

  // zero out all remaining constants
  GLOBAL_num_files = 0;
  GLOBAL_nfiles = 0;

  return;
}

static inline char get_buffchar(char *buf, size_t bufsize, size_t *cur_ind,
                                size_t *remaining, FILE *stream){
  if(*cur_ind == *remaining){
    *cur_ind = 0;
    *remaining = fread(buf, 1, bufsize, stream);
  }
  return *remaining ? buf[(*cur_ind)++] : 0;
}

static void truncate_file(const char* fname, size_t size){
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

static int nohash_name_cmpfunc(const void *a, const void *b){
  // same as above, but skip the hash comparison
  const char *aa = *(const char **)a;
  const char *bb = *(const char **)b;

  // sort first by string length
  int v1 = strlen(aa), v2 = strlen(bb);
  if (v1 != v2) return v1 - v2;

  return strcmp(aa, bb);
}

static int leaf_index_compar(const void *a, const void *b){
  l_uint aa = GLOBAL_leaf[*(l_uint*)a]->count;
  l_uint bb = GLOBAL_leaf[*(l_uint*)b]->count;
  return (aa > bb) - (aa < bb);
}

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

/**************************/
/* File Mergesort Helpers */
/**************************/

static void kway_mergesort_file_inplace(const char* f1, l_uint nlines,
                          size_t element_size,
                          l_uint block_size, int buf_size,
                          int num_bins, int output_size,
                          int (*compar)(const void *, const void *),
                          const int verbose){
  /*
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
  FILE *fileptr;

  LoserTree *mergetree = LT_alloc(num_bins, output_size, element_size, compar);
  GLOBAL_mergetree = mergetree;
  size_t cur_start, to_read;
  int empty_bin, LT_total_bins;
  l_uint nblocks, num_iter;
  double prev_progress, cur_progress;
  l_uint cur_fsize, max_fsize=0;

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
    //file2 = safe_fopen(cur_target, "wb");
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
    rewind(fileptr);
  }
  // cur_source will always be the file we just WROTE to here

  if(verbose == VERBOSE_ALL){
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

static void kway_mergesort_file(const char* f1, const char* f2, l_uint nlines,
                          size_t element_size,
                          l_uint block_size, int buf_size,
                          int num_bins, int output_size,
                          int (*compar)(const void *, const void *),
                          const int verbose){
  /*
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
  l_uint cur_fsize, max_fsize=0;

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

    fclose_tracked(2);
    tmp_swap_char = cur_source;
    cur_source = cur_target;
    cur_target = tmp_swap_char;
  }
  // cur_source will always be the file we just WROTE to here

  if(verbose == VERBOSE_ALL){
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
   * After sorting, file is a bunch of {long, double} values
   * We need to discard the long values and replace the file
   * with just the double values. nfilename will become the
   * neighbors file, and wfilename will become the weights file.
   */

  double *valuebuffer = safe_malloc(L_SIZE * FILE_READ_CACHE_SIZE);
  edge *edgebuffer = safe_malloc(sizeof(edge)*FILE_READ_CACHE_SIZE);

  FILE *f_r, *f_w;
  f_r = safe_fopen(nfilename, "rb");
  f_w = safe_fopen(nfilename, "rb+"); // wb will delete the file

  l_uint cur_line = 0;
  l_uint nread = 0;

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

static void unique_strings_with_sideeffects(char **names, int num_to_sort,
                                            int *InsertPoint, uint *counts,
                                            int useCounts){
  /*
   * This code is duplicated a lot, so just putting it here for consistency
   * This uniques the set of strings **names and stores additional information
   *  -     InsertPoint: number of unique strings
   *  -          counts: number of each unique string (ignored if !useCounts)
   */
  int insert_point = 0;
  uint cur_len;

  // first sort the array
  qsort(names, num_to_sort, sizeof(char*), nohash_name_cmpfunc);

  // next, unique the values
  if(useCounts) counts[0] = !!names[0][MAX_NODE_NAME_SIZE];
  cur_len = strlen(names[0]);
  for(int i=1; i<num_to_sort; i++){
    if(cur_len != strlen(names[i]) || (strcmp(names[i], names[insert_point]) != 0)){
      // if the string is different, save it
      insert_point++;
      if(useCounts) counts[insert_point] = 1;
      if(insert_point != i){
        memcpy(names[insert_point], names[i], MAX_NODE_NAME_SIZE+1);
        names[insert_point][MAX_NODE_NAME_SIZE] = 0;
      }
      cur_len = strlen(names[insert_point]);
    } else if(useCounts) {
      // else it's the same, so increment the corresponding count
      // (last byte stores if it's an edge that should be counted)
      counts[insert_point] += !!names[i][MAX_NODE_NAME_SIZE];
    }
  }
  insert_point++;

  // side effect writes
  *InsertPoint = insert_point;
  return;
}

static h_uint hash_file_vnames_trie(const char* fname, prefix *trie, h_uint next_index,
  const char sep, const char line_sep, int v, int is_undirected, int ignore_weights){
  /*
   * fname: .tsv list of edges
   * dname: directory of hash codes
   * hashfname: file to write vertices to
   */
  FILE *f = safe_fopen(fname, "rb");
  size_t rc_size = FILE_READ_CACHE_SIZE*8;
  size_t remaining=0, rcache_i=0;
  // size + 1 so that there's space for marking if it's an outgoing edge
  char *vname;
  char *wbuf = safe_malloc(sizeof(char) * NODE_NAME_CACHE_SIZE);
  char **restrict namecache = safe_malloc(sizeof(char*) * NODE_NAME_CACHE_SIZE);
  uint *restrict str_counts = safe_malloc(sizeof(uint) * NODE_NAME_CACHE_SIZE);
  char *read_cache = safe_malloc(rc_size);
  int num_unique;
  double weight, max_weight = ignore_weights ? 1.0 : 0.0;

  for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) namecache[i] = safe_malloc(MAX_NODE_NAME_SIZE+1);

  int cur_pos = 0, cachectr=0;
  char c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, f);
  l_uint print_counter = 0;

  if(v >= VERBOSE_BASIC) Rprintf("\tReading file %s...\n", fname);

  while(!feof(f) || remaining){
    // going to assume we're at the beginning of a line
    // lines should be of the form `start end weight` or `start end`
    // separation is by char `sep`
    for(int iter=0; iter<2; iter++){
      vname = namecache[cachectr];
      memset(vname, 0, MAX_NODE_NAME_SIZE+1);
      cur_pos = 0;
      while(c != sep && c != line_sep){
        vname[cur_pos++] = c;
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, f);
        if(cur_pos == MAX_NODE_NAME_SIZE-1) // max size has to include the null terminator
          error("Node name is larger than max allowed name size.\n");

        if(feof(f) && !remaining) error("Unexpected end of file.\n");
      }

      // mark if edge is outgoing -- always if undirected, otherwise only if the first node name
      vname[MAX_NODE_NAME_SIZE] = is_undirected ? 1 : !iter;
      str_counts[cachectr] = vname[MAX_NODE_NAME_SIZE];

      cachectr++;
      if(cachectr == NODE_NAME_CACHE_SIZE){
        unique_strings_with_sideeffects(namecache, cachectr, &num_unique, str_counts, TRUE);
        for(int i=0; i<num_unique; i++)
          next_index = insert_into_trie(namecache[i], trie, next_index, str_counts[i]);
        cachectr = 0;
      }

      // if lines are of the form `start end`, we need to leave c on the terminator
      if(c == sep)
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, f);
    }
    if(!ignore_weights){
      // read in the weight
      cur_pos = 0;
      memset(wbuf, 0, MAX_NODE_NAME_SIZE);
      while(c != line_sep && (!feof(f) || remaining)){
        wbuf[cur_pos++] = c;
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, f);
      }
      weight = atof(wbuf);
      if(weight > max_weight) max_weight = weight;
    } else {
      max_weight = 1.0;
      while(c != line_sep && (!feof(f) || remaining))
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, f);
    }

    if(c == line_sep) c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, f);
    print_counter++;
    if(!(print_counter % PRINT_COUNTER_MOD)){
      if(v == VERBOSE_ALL) Rprintf("\t%" lu_fprint " lines read\r", print_counter);
      else R_CheckUserInterrupt();
    }
  }

  if(cachectr){
    unique_strings_with_sideeffects(namecache, cachectr, &num_unique, str_counts, TRUE);
    for(int i=0; i<num_unique; i++)
        next_index = insert_into_trie(namecache[i], trie, next_index, str_counts[i]);
  }

  if(v >= VERBOSE_BASIC) Rprintf("\t%" lu_fprint " lines read\n", print_counter);
  fclose_tracked(1);
  GLOBAL_max_weight = max_weight;

  for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) free(namecache[i]);
  free(namecache);
  free(read_cache);
  free(wbuf);
  return next_index;
}

static l_uint reindex_trie_and_write_counts(prefix *trie,
                                            l_uint max_seen, int verbose,
                                            l_uint* nseen){
  // Vertices are 0-indexed
  // assume that we've already opened the file, since we'll call this recursively
  // I'm NOT going to reindex, since the read indices will be better for cache locality
  // I am going to return the max count we see, for an auto-determination for max_iterations
  if(!trie) return max_seen;
  uint8_t bits_remaining = trie->count1 + trie->count2;
  uint8_t ctr = 0;
  l_uint cur_index;
  if(trie->bmap1 & 1){
    leaf *l = (leaf*)(trie->child_nodes[ctr++]);
    cur_index = l->index;
    GLOBAL_all_leaves[cur_index] = l;
    // save the maximum degree we observe
    max_seen = max_seen > l->count ? max_seen : l->count;

    l->edge_start = l->count; // will handle this later

    // set count (cluster) to index+1
    l->count = cur_index+1;
    (*nseen)++;
    if(((*nseen)+1) % PRINT_COUNTER_MOD == 0){
      R_CheckUserInterrupt();
      if(verbose == VERBOSE_ALL) Rprintf("\tProcessed %" lu_fprint " vertices\r", *nseen);
    }
  }
  while(ctr < bits_remaining)
    max_seen = reindex_trie_and_write_counts(trie->child_nodes[ctr++],
                                              max_seen, verbose, nseen);

  return max_seen;
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

static l_uint csr_compress_edgelist_trie_batch(const char* edgefile, prefix *trie,
                                  const char* fweight, const char* fneighbor,
                                  const char sep, const char linesep, l_uint num_v, int v,
                                  const int is_undirected,
                                  const int ignore_weights, const int sort_inplace){
  /*
   * This should be called after we've already read in all our files
   * critically, ensure we're rewritten our ftable file such that it is
   * cumulative counts and not vertex counts
   *
   * Error checking can be reduced because we would have caught it earlier
   *
   * File inputs:
   *  fweight: file where edge weights will be written
   *  fneighbors: file where edge end nodes will be written
   *
   * Changes from previous version:
   *  use fread(...) to read into buffer rather than getc.
   */

  int cachectr = 0;
  char *read_cache;
  char *restrict vname;
  edge *edges;
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
  vname = safe_malloc(MAX_NODE_NAME_SIZE);
  edges = safe_malloc(FILE_READ_CACHE_SIZE * sizeof(edge));

  FILE *edgelist, *neighbortable;
  edgelist = safe_fopen(edgefile, "rb");

  // open neighbors table
  neighbortable = safe_fopen(fneighbor, "ab");

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
        vname[stringctr++] = c;
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
      }
      // short circuit to skip a length-0 name
      if(stringctr == 0){
        i--;
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
        continue;
      }
      vname[stringctr] = 0;
      inds[i] = find_index_for_prefix(vname, trie);

      // advance one past the separator if it isn't linesep
      // it would equal linesep if we don't have weights included
      if(c == sep)
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
    }
    stringctr = 0;
    if(!ignore_weights){
      // read in the weight
      memset(vname, 0, MAX_NODE_NAME_SIZE);
      c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);

      // need the double check in case user forgets a trailing \n
      while(c != linesep && (remaining || !feof(edgelist))){
        if(stringctr == MAX_NODE_NAME_SIZE)
          error("Incomplete entry read (suspect line %" lu_fprint ")", print_counter+1);
        vname[stringctr++] = c;
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
      }
      weight = atof(vname);
    } else {
      weight = 1.0;
      while(c != linesep && (remaining || !feof(edgelist)))
        c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
    }

    // add to cache
    edges[cachectr++] = compressEdgeValues(inds[0], inds[1], weight);
    if(is_undirected)
      edges[cachectr++] = compressEdgeValues(inds[1], inds[0], weight);

    // advance one past the separator
    c = get_buffchar(read_cache, rc_size, &rcache_i, &remaining, edgelist);
    if((cachectr+is_undirected) >= FILE_READ_CACHE_SIZE || (feof(edgelist) && !remaining)){
      // sort the block and write it to the neighbors file
      qsort(edges, cachectr, sizeof(edge), edge_compar);
      nedges += safe_fwrite(edges, sizeof(edge), cachectr, neighbortable);
      cachectr = 0;
    }

    print_counter++;
    if(!(print_counter % PRINT_COUNTER_MOD)){
      if(v == VERBOSE_ALL) Rprintf("\t%" lu_fprint " edges read\r", print_counter);
      R_CheckUserInterrupt();
    }
  }

  if(v >= VERBOSE_BASIC) Rprintf("\t%" lu_fprint " edges read\n", print_counter);
  fclose_tracked(2);
  free(edges);
  free(vname);
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
    if(tmp_cl == new_clust || weights[i] < CLUSTER_MIN_WEIGHT) continue;
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
                          float self_loop_weight,
                          FILE *weightsfile, FILE *neighborfile,
                          ArrayQueue *queue, float atten_param){
  /*
   * Determine number of edges using the table file (next - cur)
   * Inputs:
   *  -              ind: 0-indexed vertex id
   *  - self_loop_weight: self loop weight
   *  -          offsets: file pointer to offsets in CSR format
   *  -      weightsfile: file pointer to weights for each edge
   *  -     neighborfile: file pointer to end nodes for each edge
   *  -            queue: queue of nodes to process
   *  -      atten_param: current attenuation parameter
   *
   * Self loops are just added as an additional edge, so the total
   * number of edges is num_edges+1
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
  int weight_sign = 1;
  for(l_uint i=0; i<num_edges; i++){
    neighbors[i] = GLOBAL_all_leaves[indices[i]];
    // attenuate edges, using adaptive scaling
    // see https://doi.org/10.1103/PhysRevE.83.036103, eqns 4-5
    sufficient_weight[i] = weights_arr[i] >= self_loop_weight;

    // since we're going to allow negative weights, I'm going to make sure
    // that negatives are always negative if included
    weight_sign = weights_arr[i] < 0 ? -1 : 1;
    weights_arr[i] *= 1-(atten_param * neighbors[i]->dist);
    if(weight_sign < 0 && weights_arr[i] > 0) weights_arr[i] *= -1;
    if(fabs(weights_arr[i]) < CLUSTER_MIN_WEIGHT) weights_arr[i] = 0;
  }

  // sort both leaves and weights by assigned cluster
  // note that this just sorts the indexes vector, not the array itself
  for(l_uint i=0; i<num_edges; i++) indices[i] = i;
  GLOBAL_leaf = neighbors;
  qsort(indices, num_edges, L_SIZE, leaf_index_compar);


  // figure out the cluster to update to
  double max_weight=0, cur_weight=0, cur_error=0;
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
      min_dist = -1;
      cur_error = 0;
      found_sufficient_weight = 0;
    }
    if(weights_arr[indices[i]]){
      kahan_accu(&cur_weight, &cur_error, weights_arr[indices[i]]);
      min_dist = min_dist > cur_neighbor->dist ? cur_neighbor->dist : min_dist;
      if(sufficient_weight[indices[i]])
        found_sufficient_weight = 1;
    }
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
    // (due to attenuation), or if none are larger than self-loop
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
                          const float self_loop_weight, const double atten_pow){
  GLOBAL_verts_remaining = num_v;

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
  /*
   * TODO: we can parallelize this, right?
   *    - set up the number of threads
   *    - each thread can act on nodes, order is random but that shouldn't affect results
   *    - likely need multiple file pointers, but file is read only so that should be ok
   *    - have a shared lock that cycles through threads for the order they write in
   *    - concerns on reproducibility if traversal order is a race condition
   *    - concerns on performance benefit if the main bottleneck is disk reads
   */
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
    update_node_cluster(next_vert, self_loop_weight,
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
                            const float self_loop_weight, const double atten_pow,
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
                  num_v, num_iter, v, self_loop_weight, atten_pow);

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
  cluster_file(weightsfile, neighborfile, num_v, num_iter, v, self_loop_weight, atten_pow);

  return;
}

static l_uint write_output_clusters_trie(FILE *outfile, prefix *trie, l_uint *clust_mapping,
                                char *s, int cur_pos, char *write_buf, const size_t num_bytes,
                                const char *seps, l_uint num_v, int verbose){
  // we re-indexed the trie according to a DFS
  // thus if we just traverse the same way, we'll get the names in order
  // the hardest part here is ensuring we keep track of the character string
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

SEXP R_LPOOM_cluster(SEXP FILENAME, SEXP NUM_EFILES, // files
                    SEXP OUTDIR, SEXP OUTFILES,  // more files
                    SEXP SEPS, SEXP ITER, SEXP VERBOSE, // control flow
                    SEXP IS_UNDIRECTED, SEXP ADD_SELF_LOOPS, // optional adjustments
                    SEXP IGNORE_WEIGHTS, SEXP CONSENSUS_WEIGHTS,
                    SEXP SORT_INPLACE, SEXP ATTEN_POWER){
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
   *
   * Rough runtimes on my machine:
   *  - reading nodes: ~1,000,000 lines / sec.
   *  - reading edges: ~2,000,000 lines / sec.
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
  GLOBAL_max_weight = 1.0;

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

  // consensus stuff
  const int consensus_len = length(CONSENSUS_WEIGHTS);
  const double* consensus_w = REAL(CONSENSUS_WEIGHTS);

  // timing
  time_t time1, time2;

  // first hash all vertex names
  time1 = time(NULL);
  if(verbose >= VERBOSE_BASIC) Rprintf("Building trie for vertex names...\n");
  for(int i=0; i<num_edgefiles; i++){
    edgefile = CHAR(STRING_ELT(FILENAME, i));
    num_v = hash_file_vnames_trie(edgefile, GLOBAL_trie, num_v, seps[0], seps[1],
                                  verbose, is_undirected, ignore_weights);
  }
  time2 = time(NULL);
  if(verbose >= VERBOSE_BASIC) report_time(time1, time2, "\t");

  // allocate space for leaf counters (one extra space to hold final value of cumulative counts)
  GLOBAL_all_leaves = safe_malloc(sizeof(leaf *) * (num_v+1));

  // next, reformat the file to get final counts for each node
  if(verbose >= VERBOSE_BASIC) Rprintf("Tidying up internal tables...\n");
  l_uint print_val = 0;
  l_uint max_degree = reindex_trie_and_write_counts(GLOBAL_trie, 0, verbose, &print_val);

  // change edge_start values to cumulative counts
  GLOBAL_all_leaves[num_v] = malloc(sizeof(leaf));
  GLOBAL_all_leaves[num_v]->edge_start = 0;
  l_uint running_sum = 0;
  for(l_uint i=0; i<=num_v; i++){
    leaf* tmp_leaf = GLOBAL_all_leaves[i];
    running_sum += tmp_leaf->edge_start;
    tmp_leaf->edge_start = running_sum - tmp_leaf->edge_start;
  }

  if(verbose >= VERBOSE_BASIC) Rprintf("\tFound %" lu_fprint " unique vertices!\n", num_v);

  // get base number of iterations
  max_degree = (l_uint)(sqrt((double)max_degree));
  size_t num_bits = sizeof(aq_int) * 8 - 1; // signed, so we have one less bit to work with
  base_iter = ((aq_int)(1)) << num_bits;
  base_iter = base_iter > max_degree ? max_degree : base_iter;
  base_iter = base_iter < 5 ? 5 : base_iter; // minimum of 5 iterations per node
  int filled_iter = 0;

  for(int i=0; i<num_ofiles; i++){
    if(num_iter[i] == 0){
      num_iter[i] = base_iter;
      filled_iter = 1;
    }
  }
  if(filled_iter && verbose >= VERBOSE_BASIC)
    Rprintf("\tAutomatically setting zero-value iterations to %d\n", base_iter);

  // then, we'll create the CSR compression of all our edges
  time1 = time(NULL);
  if(verbose >= VERBOSE_BASIC) Rprintf("Reading in edges...\n");
  l_uint num_edges = 0;
  for(int i=0; i<num_edgefiles; i++){
    edgefile = CHAR(STRING_ELT(FILENAME, i));
    num_edges += csr_compress_edgelist_trie_batch(edgefile, GLOBAL_trie,
                                weightsfile, neighborfile,
                                seps[0], seps[1], num_v, verbose,
                                is_undirected,
                                ignore_weights, use_inplace_sort);
  }
  // sort the file and split it into neighbor and weight
  if(FILE_READ_CACHE_SIZE < num_edges){
    if(use_inplace_sort){
      kway_mergesort_file_inplace(neighborfile,
                                  num_edges, sizeof(edge),
                                  FILE_READ_CACHE_SIZE,
                                  MERGE_INPUT_SIZE,
                                  MAX_BINS_FOR_MERGE, // bins to merge with
                                  MERGE_OUTPUT_SIZE, // size of output buffer
                                  edge_compar, verbose);
    } else {
      kway_mergesort_file(neighborfile, weightsfile,
                                  num_edges, sizeof(edge),
                                  FILE_READ_CACHE_SIZE,
                                  MERGE_INPUT_SIZE,
                                  MAX_BINS_FOR_MERGE, // bins to merge with
                                  MERGE_OUTPUT_SIZE, // size of output buffer
                                  edge_compar, verbose);
    }
  }
  split_sorted_file(neighborfile, weightsfile, num_edges, verbose);
  time2 = time(NULL);
  //check_mergedsplit_file(neighborfile, weightsfile);

  if(verbose >= VERBOSE_BASIC) report_time(time1, time2, "\t");

  char vert_name_holder[MAX_NODE_NAME_SIZE];
  char write_buffer[PATH_MAX];
  for(int i=0; i<num_ofiles; i++){
    reset_trie_clusters(num_v);
    time1 = time(NULL);
    if(consensus_len){
      consensus_cluster_oom(weightsfile, neighborfile, dir, num_v,
                            num_iter[i], verbose, self_loop_weights[i], atten_power[i],
                            consensus_w, consensus_len);

    } else {
      if(verbose >= VERBOSE_BASIC) Rprintf("Clustering...\n");
      cluster_file(weightsfile, neighborfile, num_v, num_iter[i], verbose,
                    self_loop_weights[i], atten_power[i]);
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
  UNPROTECT(1);
  return RETURN_VAL;
}
