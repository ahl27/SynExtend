#ifndef UTIL_FILE
#define UTIL_FILE

// Various utility functions
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <R_ext/Random.h>
#include <R.h>
#include <Rdefines.h>

typedef unsigned int uint;

// SEED BEFORE CALLING ANY RANDOM FUNCTION


/*** Random numbers ***/

// random float number in range [0, 1]
double inline frand(){ return unif_rand(); }

// random integer
int inline irand(){ return (int) floor(frand() * INT_MAX);}

// random number from normal distribution with mean mu and standard deviation sd
double inline rnorm(double mu, double sd){ return sd * (norm_rand()) + mu; }


/*** Random permutations using Fisher-Yates shuffling ***/

// Returns an array with a random sample of [0,n)
int *sample(int n);

// Shuffles an array of integers
void shuffle_int_(int *x, int n);

// Shuffles an array of unsigned integers
void shuffle_uint_(uint *x, int n);

// Shuffles an array of doubles
void shuffle_double_(double *x, int n);

// Shuffles an array of char
void shuffle_char_(char *x, int n);

#define shuffle(xtype, x, n) shuffle_ ## xtype ##_(x, n)

/*** Random number generator using xorshift+ ***/

// 64 bit generation //
struct RNGstate64 {
  uint64_t state[2];
};

void seedRNGState64(struct RNGstate64 *r, uint64_t seed);

uint64_t xorshift128p(struct RNGstate64 *r);

// 32 bit generation //
struct RNGstate32 {
  uint32_t state[4];
};

static inline uint32_t rotl(const uint32_t x, int k) {
  return (x << k) | (x >> (32 - k));
}

void seedRNGState32(struct RNGstate32 *r, uint64_t seed);

uint32_t xorshift32b(struct RNGstate32 *r);

/*** Random math functions ***/
long inline doubleFactorial(int n){
  long retval = 1;
  while (n > 0) retval *= n--;
  return retval;
}

/*** Dendrapply substitutes ***/
void rdendrapplyhelper(SEXP node, SEXP f, SEXP env);
SEXP rdendrapply(SEXP tree, SEXP fn, SEXP env);
void rpdendrapplyhelper(SEXP node, SEXP f, SEXP env);
SEXP rpdendrapply(SEXP tree, SEXP fn, SEXP env);


/*** Other Functions ***/
void genCostMatrix(double *m1, double *m2, int *nc1p, int *nc2p, int *nrp, double *costMat, int *idxLookup);
void R_combineDistObj(double *d1, double *d2, int *pos1, int *pos2, int *n1, int *n2);

#endif