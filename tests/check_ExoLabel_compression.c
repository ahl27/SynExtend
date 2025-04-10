#include<stdint.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#define h_uint uint64_t
#define uint uint_fast32_t
#define strlen_uint uint_least16_t
#define l_uint uint64_t
#define lu_fprint PRIu64
#define w_float float
#define aq_int int16_t // size of "seen" counter in ArrayQueue

static const int BITS_FOR_WEIGHT = 16;
static const int BITS_AFTER_DECIMAL = 12;
static const l_uint MAX_NUM_NODES = (1ULL << (64-BITS_FOR_WEIGHT)) - 1;
static const l_uint MAX_POSSIBLE_WEIGHT = (1ULL << BITS_FOR_WEIGHT) - 1;
static const double COMPRESSED_CONVERTER = 1ULL << BITS_AFTER_DECIMAL;

typedef struct {
  l_uint v;
  l_uint w;
} edge;

static edge compressEdgeValues(l_uint v1, l_uint v2, double weight){
  /*
   * Using 4.12 compression for weights and 48-bit representation for v2
   * See line 103 for an explanation of how these were chosen.
   *
   * Note that weights are guaranteed to be positive
   */
  edge e;
  e.v = v1;
  l_uint v2_comp = v2;
  v2_comp <<= BITS_FOR_WEIGHT;

  weight = log2(weight+1);
  weight *= COMPRESSED_CONVERTER;
  l_uint w_comp = weight > MAX_POSSIBLE_WEIGHT ? MAX_POSSIBLE_WEIGHT : floor(weight);

  l_uint max_weight_bits = (1UL << BITS_FOR_WEIGHT) - 1;

  w_comp = w_comp & MAX_POSSIBLE_WEIGHT;

  e.w = v2_comp | w_comp;
  return e;
}

static void decompressEdgeValue(l_uint compressed, l_uint *v2, w_float *w){
  /*
   * See compressEdgeValues for description
   * 64-bit unsigned integer
   * first 48 bits are the index
   * last 16 bits are the compressed weight (4.12 format)
   */
  // mask to get the weight out
  l_uint max_weight_bits = (1UL << BITS_FOR_WEIGHT) - 1;

  // get the weight in compressed form
  double w_comp = (double) ((l_uint)(compressed & max_weight_bits));
  w_comp /= COMPRESSED_CONVERTER;
  w_comp = pow(2, w_comp) - 1;

  *w = (w_float)w_comp;
  *v2 = compressed >> BITS_FOR_WEIGHT;

  return;
}


int main(){
  l_uint v;
  w_float w;

  l_uint ind = MAX_NUM_NODES;
  printf("Maximum node check:\n");
  for(int i=0; i<10; i++){
    double weight = (double)i*10;
    decompressEdgeValue(compressEdgeValues(ind, ind, weight).w, &v, &w);
    if(v != ind){
      printf("%llu %llu %.5f %.5f (diff %.5f)\n", ind, v, weight, w, fabs(w-weight));
    }
  }
  printf("Done.\n\n");


  printf("Accuracy in 0-1:\n");
  for(int i=0; i<1000; i++){
    double weight = (double)i / 1000;
    decompressEdgeValue(compressEdgeValues(1, 1, weight).w, &v, &w);
    if(fabs(w-weight) > 0.0003){
      printf("%.5f %.5f (diff %.5f)\n", weight, w, fabs(w-weight));
    }
  }
  printf("Done.\n\n");

  printf("Accuracy in 60,000-65,000:\n");
  for(int i=0; i<5; i++){
    double weight = (double)i * 1000 + 60000;
    decompressEdgeValue(compressEdgeValues(1, 1, weight).w, &v, &w);
    if(fabs(w-weight) > 0.0003){
      printf("%.5f %.5f (diff %.5f)\n", weight, w, fabs(w-weight));
    }
  }
  printf("Done.\n");
}
