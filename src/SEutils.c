#include "SEutils.h"

void *safe_malloc(size_t size){
  void *data = malloc(size);
  if(!data) error("Could not allocate %zu bytes (do you have enough RAM?)", size);
  return data;
}

void *safe_calloc(size_t nitems, size_t size){
  void *data = calloc(nitems, size);
  if(!data) error("Could not allocate %zu bytes (do you have enough RAM?)", size);
  return data;
}

void *safe_realloc(void *ptr, size_t new_size){
  void *tmp = ptr;
  ptr = realloc(ptr, new_size);
  if(!ptr){
    error("Could not re-allocate %zu bytes (do you have enough RAM?)", new_size);
    if(tmp) free(tmp);
  }

  return ptr;
}

int *sample(int n){
  int *r = malloc(sizeof(int) * n);
  int j;
  for (int i=0; i<n; i++){
    j = irand() % (i+1);
    if (j != i)
      r[i] = r[j];
    r[j] = i;
  }

  return r;
}

void shuffle_int_(int *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

void shuffle_uint_(uint *x, int n){
  int j;
  uint tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}


void shuffle_double_(double *x, int n){
  int j;
  double tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

void shuffle_char_(char *x, int n){
  int j;
  char tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

void seedRNGState64(struct RNGstate64 *r, uint64_t seed){
  // splitmix64 initialization, see https://en.wikipedia.org/wiki/Xorshift
  seed += 0x9E3779B97F4A7C15;
  seed = (seed ^ (seed >> 30)) * 0xBF58476D1CE4E5B9;
  seed = (seed ^ (seed >> 27)) * 0x94D049BB133111EB;
  seed ^= (seed >> 31);
  r->state[0] = (uint32_t) seed;
  r->state[1] = (uint32_t) (seed >> 32);
  return;
}

uint64_t xorshift128p(struct RNGstate64 *r){
  uint64_t t = r->state[0];
  uint64_t const s = r->state[1];

  r->state[0] = s;
  t ^= (t << 23);
  t ^= (t >> 18);
  t ^= s ^ (s >> 5);
  r->state[1] = t;

  return t+s;
}

void seedRNGState32(struct RNGstate32 *r, uint64_t seed){
  // splitmix64 initialization, see https://en.wikipedia.org/wiki/Xorshift
  seed += 0x9E3779B97F4A7C15;
  seed = (seed ^ (seed >> 30)) * 0xBF58476D1CE4E5B9;
  seed = (seed ^ (seed >> 27)) * 0x94D049BB133111EB;
  seed ^= (seed >> 31);
  r->state[0] = (uint32_t) seed;
  r->state[1] = (uint32_t) (seed >> 32);

  seed += 0x9E3779B97F4A7C15;
  seed = (seed ^ (seed >> 30)) * 0xBF58476D1CE4E5B9;
  seed = (seed ^ (seed >> 27)) * 0x94D049BB133111EB;
  seed ^= (seed >> 31);
  r->state[2] = (uint32_t) seed;
  r->state[3] = (uint32_t) (seed >> 32);
  return;
}

// adapted from https://prng.di.unimi.it/xoshiro128plus.c
uint32_t xorshift32b(struct RNGstate32 *rngstate) {
  uint32_t *r = rngstate->state;
  const uint32_t result = r[0] + r[3];
  const uint32_t t = r[1] << 9;

  r[2] ^= r[0];
  r[3] ^= r[1];
  r[1] ^= r[2];
  r[0] ^= r[3];

  r[2] ^= t;

  r[3] = rotl(r[3], 11);

  return result;
}


/* needs C11 :/
#define shuffle(x, n) _Generic((x),
  int *: shuffle_i(x, n),
  double *: shuffle_d(x, n),
  char *: shuffle_c(x, n))
*/

/* Dendrapply Methods */

// rdendrapply, should work the same as base dendrapply but faster
// possibly a deeper recursion stack as well
void rdendrapplyhelper(SEXP node, SEXP f, SEXP env){
  if (isNull(getAttrib(node, install("leaf")))){
    int n = length(node);
    for(int i=0; i<n; i++){
      SEXP call = PROTECT(LCONS(f, LCONS(VECTOR_ELT(node, i), R_NilValue)));
      SET_VECTOR_ELT(node, i, R_forceAndCall(call, 1, env));
      UNPROTECT(1);
      rdendrapplyhelper(VECTOR_ELT(node, i), f, env);
    }
  }
  return;
}

SEXP rdendrapply(SEXP tree, SEXP fn, SEXP env){
  SEXP treecopy = PROTECT(duplicate(tree));
  SEXP call = PROTECT(LCONS(fn, LCONS(treecopy, R_NilValue)));
  treecopy = PROTECT(R_forceAndCall(call, 1, env));
  rdendrapplyhelper(treecopy, fn, env);
  UNPROTECT(3);
  return treecopy;
}

// rpdendrapply, recursively appiles a function passing the parent node as an arg
void rpdendrapplyhelper(SEXP node, SEXP f, SEXP env){
  if (isNull(getAttrib(node, install("leaf")))){
    int n = length(node);
    for(int i=0; i<n; i++){
      SEXP call = PROTECT(LCONS(f, LCONS(VECTOR_ELT(node, i), LCONS(node, R_NilValue))));
      SET_VECTOR_ELT(node, i, R_forceAndCall(call, 2, env));
      UNPROTECT(1);
      rpdendrapplyhelper(VECTOR_ELT(node, i), f, env);
    }
  }
  return;
}

SEXP rpdendrapply(SEXP tree, SEXP fn, SEXP env){
  SEXP treecopy = PROTECT(duplicate(tree));
  rpdendrapplyhelper(treecopy, fn, env);
  UNPROTECT(1);
  return treecopy;
}

/* Other Random Stuff */
void genCostMatrix(double *m1, double *m2, int *nc1p, int *nc2p, int *nrp, double *costMat, int *idxLookup){
  const int nc1 = *nc1p;
  const int nc2 = *nc2p;
  const int nr = *nrp;

  double minv, tmp;
  int col1, col2, row=0;
  // remember that R matrices are column-major
  for(int i=0; i<nc1; i++){
    col1 = i*nr;
    for(int j=0; j<nc2; j++){
      col2 = j*nr;
      minv = -1;
      for (int k=0; k<nr; k++){
        tmp = m1[col1+k] + m2[col2+k];
        if(tmp < minv || minv < 0){
          minv = tmp;
          row = k+1;
        }
      }
      costMat[j*nc1+i] = minv;
      idxLookup[j*nc1+i] = row;
    }
  }

  return;
}

void R_combineDistObj(double *d1, double *d2, int *pos2, int *n1, int *n2, double *mult){
  // d2 is the shorter distance matrix that will be added into d1
  const int len1 = n1[0];
  const int len2 = n2[0];
  const int iterlen = (len2*(len2-1)) / 2;

  // remember formatting--dist objects are just vectors stored with a size attribute
  // thus col1 is d1[1:(n-1)], col2 is d1[(n+1):(n+(n-1)-1)], etc.
  // diagonal is always zero and is thus not stored

  // current row/col INDEX in second dist
  int col2 = 0;
  int row2 = 1;
  double v, w;
  int m,M;
  for(int i=0; i<iterlen; i++){
    v = d2[i];
    if(pos2[row2] < pos2[col2]){
      m = pos2[row2];
      M = pos2[col2];
    } else {
      m = pos2[col2];
      M = pos2[row2];
    }
    w = mult[row2] * mult[col2];
    v *= w;

    // If the values are the same, just skip
    if(m!=M)
      d1[len1*(m - 1) - m*(m - 1)/2 + M - m - 1] += v;

    row2++;
    if(row2 == len2){
      col2++;
      row2=col2+1;
    }
  }

  return;
}




