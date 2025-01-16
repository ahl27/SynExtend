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
#define l_uint uint_fast64_t
#define lu_fprint PRIuFAST64
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
static const int MAX_READ_RETRIES = 10;
static const int MAX_WRITE_RETRIES = 10;
static const char CONSENSUS_CSRCOPY1[] = "tmpcsr1";
static const char CONSENSUS_CLUSTER[] = "tmpclust";
static const int BITS_FOR_WEIGHT = 10;

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

static void report_time(clock_t start, clock_t end, const char* prefix){
  double elapsed_time = ((double)(end - start)) / CLOCKS_PER_SEC;
  double secs;
  int mins, hours, days;
  secs = fmod(elapsed_time, 60);

  int elapsed_secs = (int)(elapsed_time - secs);
  days = elapsed_secs / 86400;
  elapsed_secs %= 86400;

  hours = elapsed_secs / 3600;
  elapsed_secs %= 3600;

  mins = elapsed_secs / 60;

  Rprintf("%sTime difference of ", prefix);
  if(days) Rprintf("%d days, ", days);
  if(hours) Rprintf("%d hrs, ", hours);
  if(mins) Rprintf("%d mins, ", mins);
  if(mins) Rprintf("%d secs\n", (int)secs);
  else Rprintf("%0.2f secs\n", secs);
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

inline float sigmoid_transform(const float w, const double slope){
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

/**************************/
/* File Function Wrappers */
/**************************/

static size_t safe_fread(void *buffer, size_t size, size_t count, FILE *stream){
  size_t found_values = fread(buffer, size, count, stream);
  if(found_values != count){
    // two scenarios:

    // 1. read past the end of the file (throw error and return)
    if(feof(stream))
      error("%s", "Internal error: fread tried to read past end of file.\n");

    // 2. some undefined reading error (retry a few times and then return)
    for(int i=0; i<MAX_READ_RETRIES; i++){
      // if we read a partial value, reset the counter back
      if(found_values) fseek(stream, -1*((int)found_values), SEEK_CUR);

      // try to read again
      found_values = fread(buffer, size, count, stream);
      if(found_values == count) return found_values;
    }

    // otherwise throw an error
    error("Internal error: fread read %zu values (expected %zu).\n", found_values, count);
  }
  return found_values;
}

static size_t safe_fwrite(void *buffer, size_t size, size_t count, FILE *stream){
  size_t written_values = fwrite(buffer, size, count, stream);
  if(written_values != count){
    // same scenarios as in safe_fread
    if(feof(stream))
      error("%s", "Internal error: fread tried to read past end of file.\n");

    for(int i=0; i<MAX_WRITE_RETRIES; i++){
      // if we read a partial value, reset the counter back
      if(written_values) fseek(stream, -1*((int)written_values), SEEK_CUR);

      // try to read again
      written_values = fwrite(buffer, size, count, stream);
      if(written_values == count) return count;
    }

    // otherwise throw an error
    error("Internal error: fwrite wrote %zu values (expected %zu).\n", written_values, count);
  }
  return written_values;
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
  if(GLOBAL_ftracker) free(GLOBAL_ftracker);

  // free array queue
  if(GLOBAL_queue) free_array_queue(GLOBAL_queue);

  // delete any files made during runtime
  // (note: does not include outfile)
  for(int i=0; i<GLOBAL_nfiles; i++){
    remove(GLOBAL_filenames[i]);
    free(GLOBAL_filenames[i]);
  }
  if(GLOBAL_filenames) free(GLOBAL_filenames);

  if(GLOBAL_all_leaves) free(GLOBAL_all_leaves);
  free_trie(GLOBAL_trie);

  if(GLOBAL_mergebuffers){
    for(int i=0; i<GLOBAL_nbuffers; i++){
      if(GLOBAL_mergebuffers[i]) free(GLOBAL_mergebuffers[i]);
    }
    if(GLOBAL_mergebuffers) free(GLOBAL_mergebuffers);
  }

  if(GLOBAL_mergetree) LT_free(GLOBAL_mergetree);

  // NULL out all remaining pointers
  GLOBAL_num_files = 0;
  GLOBAL_nfiles = 0;
  GLOBAL_trie = NULL;
  GLOBAL_ftracker = NULL;
  GLOBAL_all_leaves = NULL;
  GLOBAL_filenames = NULL;

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
    int filehandler = _open(fname, _O_BINARY | _O_RDWR);
    retval = _chsize_s(filehandler, size);
    _close(filehandler);
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
  if(verbose) Rprintf("\tSorting with %d-way merge...\n", num_bins);

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
    if(verbose){
      Rprintf("\tIteration %llu of %llu (%6.01f%% complete, used ",
                              tmpniter, nmax_iterations, cur_progress);
      report_filesize(max_fsize);
      Rprintf(")  \r");
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
          if(verbose){
            if(tmpniter == 1){
              cur_fsize = ftell(fileptr);
              if(cur_fsize > max_fsize) max_fsize = cur_fsize;
            }
            Rprintf("\tIteration %llu of %llu (%6.01f%% complete, used ",
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
          if(verbose){
            if(tmpniter == 1){
              cur_fsize = ftell(fileptr);
              if(cur_fsize > max_fsize) max_fsize = cur_fsize;
            }
            Rprintf("\tIteration %llu of %llu (%6.01f%% complete, used ",
                              tmpniter, nmax_iterations, cur_progress);
            report_filesize(max_fsize);
            Rprintf(")  \r");
          }
        R_CheckUserInterrupt();
      }
    }

    mergetree->nwritten = 0;
    block_size *= num_bins;
    rewind(fileptr);
  }
  // cur_source will always be the file we just WROTE to here

  if(verbose) Rprintf("\n");
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
  if(verbose) Rprintf("\tSorting with %d-way merge...\n", num_bins);

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
    if(verbose){
      Rprintf("\tIteration %llu of %llu (%6.01f%% complete, used ",
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
          if(verbose){
            if(tmpniter == 1){
              cur_fsize = ftell(file1) + ftell(file2);
              if(cur_fsize > max_fsize) max_fsize = cur_fsize;
            }
            Rprintf("\tIteration %llu of %llu (%6.01f%% complete, used ",
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
        if(verbose){
          if(tmpniter == 1){
            cur_fsize = ftell(file1) + ftell(file2);
            if(cur_fsize > max_fsize) max_fsize = cur_fsize;
          }
          Rprintf("\tIteration %llu of %llu (%6.01f%% complete, used ",
                            tmpniter, nmax_iterations, cur_progress);
          report_filesize(max_fsize);
          Rprintf(")  \r");
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

  if(verbose){
    Rprintf("\tIteration %llu of %llu (%6.01f%% complete, used ",
                      tmpniter, nmax_iterations, cur_progress);
    report_filesize(max_fsize);
    Rprintf(")  \r");
  }
  R_CheckUserInterrupt();
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
    for(int i=0; i<nread; i++){
      decompressEdgeValue(valuebuffer[i], &(indices[i]), &(weights[i]));
    }

    safe_fwrite(indices, L_SIZE, nread, f_w);
    safe_fwrite(weights, W_SIZE, nread, wf_w);
    cur_line += nread;
    cur_progress = 100*((double)cur_line) / nlines;
    if(cur_progress - prev_progress > 0.5 || cur_progress == 100){
      if(v) Rprintf("\tTidying up edges (%6.01f%% Complete)\r", cur_progress);
      prev_progress = cur_progress;
      R_CheckUserInterrupt();
    }
  }
  if(v) Rprintf("\n");
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

  if(v) Rprintf("\tReading file %s...\n", fname);

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
      if(v) Rprintf("\t%" lu_fprint " lines read\r", print_counter);
      else R_CheckUserInterrupt();
    }
  }

  if(cachectr){
    unique_strings_with_sideeffects(namecache, cachectr, &num_unique, str_counts, TRUE);
    for(int i=0; i<num_unique; i++)
        next_index = insert_into_trie(namecache[i], trie, next_index, str_counts[i]);
  }

  if(v) Rprintf("\t%" lu_fprint " lines read\n", print_counter);
  fclose_tracked(1);
  GLOBAL_max_weight = max_weight;

  for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) free(namecache[i]);
  free(namecache);
  free(read_cache);
  free(wbuf);
  return next_index;
}

static l_uint reindex_trie_and_write_counts(prefix *trie, FILE* csrfile,
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

    fseek(csrfile, cur_index*L_SIZE, SEEK_SET);
    safe_fwrite(&(l->count), L_SIZE, 1, csrfile);

    // set count (cluster) to index+1
    l->count = cur_index+1;
    (*nseen)++;
    if(((*nseen)+1) % PRINT_COUNTER_MOD == 0){
      R_CheckUserInterrupt();
      if(verbose) Rprintf("\tProcessed %" lu_fprint " vertices\r", *nseen);
    }
  }
  while(ctr < bits_remaining)
    max_seen = reindex_trie_and_write_counts(trie->child_nodes[ctr++], csrfile,
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
  }

  return;
}

static void reformat_counts(const char* curcounts, const char* mastertable,
                            l_uint n_vert, int add_self_loops){
  /*
   * Creates a new table with cumulative counts
   * leaves the old table unchanged, this will act as a temporary counts file later
   */
  const uint l_size = L_SIZE;
  l_uint cumul_total = 0, curcount;
  FILE *tmptab = safe_fopen(curcounts, "rb");
  FILE *mtab = safe_fopen(mastertable, "wb+");
  //int self_loop = add_self_loops ? 1 : 0;
  for(l_uint i=0; i<n_vert; i++){

    safe_fwrite(&cumul_total, l_size, 1, mtab);
    safe_fread(&curcount, l_size, 1, tmptab);
    cumul_total += curcount;

    // also write this information to the leaves to optimize writing edges later
    GLOBAL_all_leaves[i]->count = cumul_total; // end of block for this vertex
  }

  // ending position of file
  safe_fwrite(&cumul_total, l_size, 1, mtab);
  fclose_tracked(2);
  return;
}

static void normalize_csr_edgecounts_batch(const char* ftable, const char* fweights,
                                            l_uint num_v, int verbose,
                                            float self_weight){
  float normalizer;
  l_uint start, end, num_entries;
  FILE *offsets = safe_fopen(ftable, "rb");
  FILE *weights = safe_fopen(fweights, "rb+");

  if(verbose) Rprintf("\tNodes remaining: %" lu_fprint "", num_v);

  // switching to just reading all edges in
  w_float *weights_arr = NULL;

  safe_fread(&end, L_SIZE, 1, offsets);
  for(l_uint i=0; i<num_v; i++){
    normalizer = self_weight;
    start = end;
    safe_fread(&end, L_SIZE, 1, offsets);
    num_entries = end - start;
    if(!num_entries) continue;

    // read in all the weights
    // note realloc on NULL is the same as malloc
    weights_arr = safe_realloc(weights_arr, W_SIZE*num_entries);
    safe_fread(weights_arr, W_SIZE, num_entries, weights);

    fseek(weights, (-1*num_entries) * W_SIZE, SEEK_CUR);

    // normalize the weights
    for(l_uint j=0; j<num_entries; j++)
      normalizer += fabsf(weights_arr[j]);
    normalizer = normalizer == 0 ? 1 : normalizer;
    for(l_uint j=0; j<num_entries; j++)
      weights_arr[j] /= normalizer;
    safe_fwrite(weights_arr, W_SIZE, num_entries, weights);

    if(i % PROGRESS_COUNTER_MOD == 0){
      if(verbose)
        Rprintf("\r\tNodes remaining: %" lu_fprint "                   ", num_v-i-1);
      else
        R_CheckUserInterrupt();
    }
  }

  if(verbose) Rprintf("\r\tNodes remaining: Done!               \n");

  if(weights_arr) free(weights_arr);
  fclose_tracked(2);
  return;
}

static l_uint csr_compress_edgelist_trie_batch(const char* edgefile, prefix *trie,
                                  const char* ftable, const char* fweight, const char* fneighbor,
                                  const char sep, const char linesep, l_uint num_v, int v,
                                  const int is_undirected, int has_self_loops,
                                  const int ignore_weights, const int sort_inplace){
  /*
   * This should be called after we've already read in all our files
   * critically, ensure we're rewritten our ftable file such that it is
   * cumulative counts and not vertex counts
   *
   * Error checking can be reduced because we would have caught it earlier
   *
   * File inputs:
   *  ftable: offsets file
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

  FILE *offsetstable, *edgelist, *neighbortable;
  offsetstable = safe_fopen(ftable, "rb+");
  if(!offsetstable){
    fclose_tracked(1); // handling for NULL is done in this function
    offsetstable = safe_fopen(ftable, "ab+");
    fclose_tracked(1);
    offsetstable = safe_fopen(ftable, "rb+");
  }
  edgelist = safe_fopen(edgefile, "rb");

  // open neighbors table
  neighbortable = safe_fopen(fneighbor, "ab");

  // note: weights/neighbors don't need to be initialized
  // see https://stackoverflow.com/questions/31642389/fseek-a-newly-created-file

  if(v) Rprintf("\tReading file %s...\n", edgefile);

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
      if(v) Rprintf("\t%" lu_fprint " edges read\r", print_counter);
      R_CheckUserInterrupt();
    }
  }

  if(v) Rprintf("\t%" lu_fprint " edges read\n", print_counter);
  fclose_tracked(3);
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
    for(int j=0; j<ctr; j++){
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

static void update_node_cluster(l_uint ind, double inflation, float self_loop_weight,
                          aq_int times_seen,
                          FILE *offsets, FILE *weightsfile, FILE *neighborfile,
                          ArrayQueue *queue){
  /*
   * Determine number of edges using the table file (next - cur)
   * If number of edges too large, use some sort of hash to bin edges, then rerun with less
   * Inputs:
   *  -          ind: 0-indexed vertex id
   *  -       offset: location where edges begin in mastertab
   *  -      offsets: file pointer to offsets in CSR format
   *  -  clusterings: file of current cluster numbers (0=unassigned)
   *  -  weightsfile: file pointer to weights for each edge
   *  - neighborfile: file pointer to end nodes for each edge
   *
   *  Slowest parts of this function are reading from disk
   *    - first two freads: 1-3 clocks
   *    - second two freads: 3-5 clocks
   *    - total function time: 9-13 clocks
   *  those four lines are 50-100% of the runtime
   */
  R_CheckUserInterrupt();
  l_uint start, end, num_edges;
  w_float *weights_arr;
  l_uint *clusts;
  l_uint *indices;
  leaf **neighbors;

  // the log scaling slope is a nice curve that plateaus so we don't grow to infinity
  double infl_pow = inflation;
  if(times_seen-1 > 0)
    infl_pow = pow(infl_pow, 1+log2((double)(times_seen-1)));
  self_loop_weight = pow(self_loop_weight, infl_pow);
  // move to information for the vertex and read in number of edges
  // TODO: can we put this in the trie?
  //        would have to have an extra value so we could get start/end
  //        from leaves[ind] and leaves[ind+1]
  fseek(offsets, L_SIZE*ind, SEEK_SET);
  safe_fread(&start, L_SIZE, 1, offsets);
  safe_fread(&end, L_SIZE, 1, offsets);

  num_edges = end - start;
  // if it has no edges we can't do anything
  if(!num_edges) return;

  weights_arr = safe_malloc(W_SIZE*num_edges);
  clusts = safe_malloc(L_SIZE*num_edges);
  // keep an extra space here so that we can put the current node too
  neighbors = safe_malloc(sizeof(leaf*)*(num_edges+1));

  // read in the edges
  fseek(weightsfile, start*W_SIZE, SEEK_SET);
  fseek(neighborfile, start*L_SIZE, SEEK_SET);
  safe_fread(weights_arr, W_SIZE, num_edges, weightsfile);

  // this is the offsets to read to in the cluster file
  safe_fread(clusts, L_SIZE, num_edges, neighborfile);

  // read in the clusters
  for(l_uint i=0; i<num_edges; i++){
    neighbors[i] = GLOBAL_all_leaves[clusts[i]];
  }
  free(clusts);

  // read in the current node too
  neighbors[num_edges] = GLOBAL_all_leaves[ind];

  // update the weights according to current inflation
  // we HAVE to renormalize because of thresholding
  // (meaning collapsing weights below CLUSTER_MIN_WEIGHT to 0)
  double total=self_loop_weight, err=0;
  for(l_uint i=0; i<num_edges; i++){
    weights_arr[i] = pow(weights_arr[i], infl_pow);
    kahan_accu(&total, &err, weights_arr[i]);
  }
  for(l_uint i=0; i<num_edges; i++){
    weights_arr[i] /= total;
    if(weights_arr[i] < CLUSTER_MIN_WEIGHT) weights_arr[i] = 0;
  }
  self_loop_weight /= total;
  if(self_loop_weight < CLUSTER_MIN_WEIGHT) self_loop_weight = 0;

  indices = safe_malloc(L_SIZE*num_edges);
  for(l_uint i=0; i<num_edges; i++) indices[i] = i;
  GLOBAL_leaf = neighbors;
  qsort(indices, num_edges, L_SIZE, leaf_index_compar);
  double max_weight, cur_weight=0, cur_error=0;
  l_uint max_clust=0, cur_clust=neighbors[num_edges]->count;

  if(num_edges && neighbors[indices[0]]->count == cur_clust){
    cur_weight += self_loop_weight;
    self_loop_weight = 0;
  }

  for(l_uint i=0; i<num_edges; i++){
    if(neighbors[indices[i]]->count != cur_clust){
      if(max_weight < cur_weight){
        max_weight = cur_weight;
        max_clust = cur_clust;
      }
      cur_clust = neighbors[indices[i]]->count;
      cur_weight = 0;
      if(self_loop_weight && cur_clust == neighbors[num_edges]->count){
        cur_weight = self_loop_weight;
        self_loop_weight = 0;
      }
      cur_error=0;
    }
    if(weights_arr[indices[i]])
      kahan_accu(&cur_weight, &cur_error, weights_arr[indices[i]]);
  }

  if(max_weight < cur_weight || max_clust == 0){
    // max_clust == 0 case is for when the node has no edges,
    // which would skip the above loop and assign the node to 0 (invalid)
    max_weight = cur_weight;
    max_clust = cur_clust;
  }

  // have to actually write the new cluster
  neighbors[num_edges]->count = max_clust;

  add_remaining_to_queue(max_clust, neighbors, weights_arr, num_edges, queue);

  free(weights_arr);
  free(neighbors);
  free(indices);
  return;
}

static void cluster_file(const char* offsets_fname, const char* weights_fname,
                          const char* neighbor_fname,
                          l_uint num_v, int max_iterations, int v,
                          double inflation, float self_loop_weight){
  GLOBAL_verts_remaining = num_v;

  // main runner function to cluster nodes
  FILE *offsetsfile = safe_fopen(offsets_fname, "rb");
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
  float pct_complete = 0, prev_pct=0;
  // randomly initialize queue and ctr file

  if(v) Rprintf("\tInitializing queues...");
  GLOBAL_queue = init_array_queue(num_v, max_iterations);
  if(v) Rprintf("done.\n\tClustering network:\n");

  if(v) Rprintf("\t0.0%% complete %s", progress[++statusctr%progbarlen]);
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

    print_counter++;
    if(!(print_counter % PROGRESS_COUNTER_MOD)){
      if(v){
        pct_complete = ((float)(num_v - GLOBAL_queue->length)) / num_v;
        if(pct_complete != 1 && fabs(pct_complete - prev_pct) < MIN_UPDATE_PCT){
          pct_complete = prev_pct;
        } else {
          prev_pct = pct_complete;
        }
        Rprintf("\r\t%0.1f%% complete %s", (pct_complete)*100, progress[++statusctr%progbarlen]);
      } else{
        R_CheckUserInterrupt();
      }
    }
    l_uint next_vert = array_queue_pop(GLOBAL_queue);
    // will either be negative (because just removed from queue)
    // or zero (because seen max_iteration times)
    aq_int num_seen = -1*GLOBAL_queue->seen[next_vert];
    if(!num_seen) num_seen = max_iterations;
    update_node_cluster(next_vert, inflation, self_loop_weight, num_seen,
                        offsetsfile, weightsfile, neighborfile, GLOBAL_queue);
  }
  if(v){
    if(max_iterations > 0)
      Rprintf("\r\t100%% complete!                \n");
    else
      Rprintf("\r\tComplete! (%d total iterations)     \n", i+1);
  }
  free_array_queue(GLOBAL_queue);
  GLOBAL_queue = NULL;
  fclose_tracked(3);

  return;
}

static void resolve_cluster_consensus(const char* clusters, const char* csrheader,
                                      const char* weightsfile, const char* neighborfile,
                                      l_uint num_v, int num_runs){
  // overwrite the files with new info
  FILE *clusts = safe_fopen(clusters, "rb");
  FILE *offsets = safe_fopen(csrheader, "wb+");

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
    tmp += counts[i];
    tmp_space[i] = tmp;
  }
  total_edges = tmp;

  // write cumulative counts to file
  l_uint cur_clust, tmp_ind, num_neighbors;
  tmp = 0;
  safe_fwrite(&tmp, L_SIZE, 1, offsets);
  safe_fwrite(tmp_space, L_SIZE, num_v, offsets);

  // write the edges
  rewind(clusts);
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

        fseek(offsets, cur_clust*L_SIZE, SEEK_SET);
        safe_fread(&tmp_ind, L_SIZE, 1, offsets);

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
  fseek(offsets, num_v*L_SIZE, SEEK_SET);
  safe_fread(&tmp_ind, L_SIZE, 1, offsets); // total number of edges
  while(total_edges){
    tmp = FILE_READ_CACHE_SIZE > total_edges ? total_edges : FILE_READ_CACHE_SIZE;
    safe_fwrite(w, W_SIZE, tmp, weights);
    total_edges -= tmp;
  }
  free(w);
  fclose_tracked(3);

  return;
}

static void consensus_cluster_oom(const char* csrfile, const char* weightsfile,
                            const char* neighborfile, const char* dir,
                            l_uint num_v, int num_iter, int v,
                            double inflation, float self_loop_weight,
                            const double* consensus_weights, const int consensus_len){

  /*
   * Inputs:
   *  - csrfile: file of offsets (read-only)
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
  l_uint num_edges = 0;
  FILE *csr = safe_fopen(csrfile, "rb");
  fseek(csr, -1*L_SIZE, SEEK_END);
  safe_fread(&num_edges, L_SIZE, 1, csr);
  fclose_tracked(1);

  l_uint *clusters = safe_malloc(L_SIZE * num_v);
  leaf *tmpleaf;

  dummyclust = safe_fopen(tmpclusterfile, "wb");
  // now we run clustering over consensus_len times
  for(int i=0; i<consensus_len; i++){
    if(v) Rprintf("Iteration %d of %d:\n", i+1, consensus_len);

    // modify weights according to sigmoid transformation
    if(v) Rprintf("\tTransforming edge weights...\n");
    copy_weightsfile_sig(transformedweights, weightsfile, num_edges, consensus_weights[i]);

    // reset cluster values
    reset_trie_clusters(num_v);

    // cluster with transformed weights
    cluster_file(csrfile, transformedweights, neighborfile,
                  num_v, num_iter, v, inflation, self_loop_weight);

    if(v) Rprintf("\tRecording results...\n");
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
  if(v) Rprintf("Reconciling runs...\n");
  resolve_cluster_consensus(tmpclusterfile, csrfile, weightsfile, neighborfile,
                            num_v, consensus_len);

  if(v) Rprintf("Clustering on consensus data...\n");
  reset_trie_clusters(num_v);
  cluster_file(csrfile, weightsfile, neighborfile, num_v, num_iter, v, 1.0, self_loop_weight);

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
        //fseek(clusterfile, tmpval*L_SIZE, SEEK_SET);
        //safe_fread(&tmpval, L_SIZE, 1, clusterfile);

        // prepare data for output
        s[cur_pos] = 0;
        snprintf(write_buf, num_bytes, "%s%c%" lu_fprint "%c", s, seps[0], tmpval, seps[1]);
        safe_fwrite(write_buf, 1, strlen(write_buf), outfile);
        num_v--;
        if(num_v % PROGRESS_COUNTER_MOD == 0){
          if(verbose)
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
                    SEXP OUTDIR, SEXP OUTFILE,  // more files
                    SEXP SEPS, SEXP ITER, SEXP VERBOSE, // control flow
                    SEXP IS_UNDIRECTED, SEXP ADD_SELF_LOOPS, // optional adjustments
                    SEXP IGNORE_WEIGHTS, SEXP NORMALIZE_WEIGHTS,
                    SEXP CONSENSUS_WEIGHTS, SEXP INFLATION_POW,
                    SEXP SORT_INPLACE){
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
  const char* tabfile = create_filename(dir, "tabfile.bin");
  const char* temptabfile = create_filename(dir, "temptabfile.bin");
  const char* weightsfile = create_filename(dir, "weights.bin");
  const char* neighborfile = create_filename(dir, "neighbors.bin");
  const char* outfile = CHAR(STRING_ELT(OUTFILE, 0));
  const char* edgefile;

  // required parameters
  const char* seps = CHAR(STRING_ELT(SEPS, 0));
  const int num_edgefiles = INTEGER(NUM_EFILES)[0];
  aq_int num_iter = INTEGER(ITER)[0];
  const int verbose = LOGICAL(VERBOSE)[0];
  l_uint num_v = 0;

  // optional parameters
  const int is_undirected = LOGICAL(IS_UNDIRECTED)[0];
  const float self_loop_weight = REAL(ADD_SELF_LOOPS)[0];
  const double inflation = REAL(INFLATION_POW)[0];
  const int add_self_loops = self_loop_weight > 0;
  const int should_normalize = LOGICAL(NORMALIZE_WEIGHTS)[0];
  const int ignore_weights = LOGICAL(IGNORE_WEIGHTS)[0];
  const int use_inplace_sort = LOGICAL(SORT_INPLACE)[0];

  // consensus stuff
  const int consensus_len = length(CONSENSUS_WEIGHTS);
  const double* consensus_w = REAL(CONSENSUS_WEIGHTS);

  // timing
  clock_t time1, time2;

  // first hash all vertex names
  time1 = clock();
  if(verbose) Rprintf("Building trie for vertex names...\n");
  for(int i=0; i<num_edgefiles; i++){
    edgefile = CHAR(STRING_ELT(FILENAME, i));
    num_v = hash_file_vnames_trie(edgefile, GLOBAL_trie, num_v, seps[0], seps[1],
                                  verbose, is_undirected, ignore_weights);
  }
  time2 = clock();
  if(verbose) report_time(time1, time2, "\t");

  // allocate space for leaf counters
  GLOBAL_all_leaves = safe_malloc(sizeof(leaf *) * (num_v+1));

  // next, reformat the file to get final counts for each node
  if(verbose) Rprintf("Tidying up internal tables...\n");
  //num_v = node_vertex_file_cleanup(dir, hashfile, temptabfile, hashindex, verbose);
  FILE *tempcounts = safe_fopen(temptabfile, "wb");
  if(!tempcounts) error("could not open file %s\n", temptabfile);

  l_uint print_val = 0;
  l_uint max_degree = reindex_trie_and_write_counts(GLOBAL_trie, tempcounts, 0, verbose, &print_val);
  fclose_tracked(1);

  if(verbose) Rprintf("\tFound %" lu_fprint " unique vertices!\n", num_v);
  if(!num_iter){
    max_degree = (l_uint)(sqrt((double)max_degree)) + 1 + add_self_loops;
    size_t num_bits = sizeof(aq_int) * 8 - 1; // signed, so we have one less bit to work with
    num_iter = ((aq_int)(1)) << num_bits;
    num_iter = num_iter > max_degree ? max_degree : num_iter;
    if(verbose) Rprintf("\tAutomatically setting iterations to %d\n", num_iter);
  }

  // next, create an indexed table file of where data for each vertex will be located
  if(verbose) Rprintf("Reformatting vertex degree file...\n");
  reformat_counts(temptabfile, tabfile, num_v, add_self_loops);

  // then, we'll create the CSR compression of all our edges
  time1 = clock();
  if(verbose) Rprintf("Reading in edges...\n");
  l_uint num_edges = 0;
  for(int i=0; i<num_edgefiles; i++){
    edgefile = CHAR(STRING_ELT(FILENAME, i));
    num_edges += csr_compress_edgelist_trie_batch(edgefile, GLOBAL_trie,
                                tabfile, weightsfile, neighborfile,
                                seps[0], seps[1], num_v, verbose,
                                is_undirected, add_self_loops,
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
  time2 = clock();
  //check_mergedsplit_file(neighborfile, weightsfile);

  if(verbose) report_time(time1, time2, "\t");

  if(verbose && should_normalize)
    Rprintf("Normalizing edge weights...\n");

  if(should_normalize)
    normalize_csr_edgecounts_batch(tabfile, weightsfile, num_v, verbose, self_loop_weight);

  reset_trie_clusters(num_v);

  time1 = clock();
  if(consensus_len){
    consensus_cluster_oom(tabfile, weightsfile, neighborfile, dir, num_v,
                          num_iter, verbose, inflation, self_loop_weight,
                          consensus_w, consensus_len);

  } else {
    if(verbose) Rprintf("Clustering...\n");
    cluster_file(tabfile, weightsfile, neighborfile, num_v, num_iter, verbose,
                  inflation, self_loop_weight);
  }
  time2 = clock();
  if(verbose) report_time(time1, time2, "\t");

  // remove these files here so that we free up the space
  remove(tabfile);
  remove(weightsfile);
  remove(neighborfile);

  // have to allocate resources for writing out
  FILE *results = safe_fopen(outfile, "wb");
  if(!results) error("Failed to open output file.");
  char vert_name_holder[MAX_NODE_NAME_SIZE];
  char write_buffer[PATH_MAX];
  l_uint *clust_mapping = safe_calloc(num_v, L_SIZE);
  CLUST_MAP_CTR = 1;
  if(verbose) Rprintf("Writing clusters to file...\n\tVertices remaining: %" lu_fprint "", num_v);
  write_output_clusters_trie(results, GLOBAL_trie, clust_mapping, vert_name_holder,
                              0, write_buffer, (size_t)PATH_MAX, seps, num_v, verbose);
  if(verbose) Rprintf("\r\tNodes remaining: Done!               \n");
  free(clust_mapping);
  fclose_tracked(1);

  cleanup_ondisklp_global_values();
  return R_NilValue;
}
