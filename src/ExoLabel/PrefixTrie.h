#ifndef PREFIXTRIE_H
#define PREFIXTRIE_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "../SEutils.h"

// testing includes
//#include <stdio.h>

#define trie_uint uint64_t
#define dist_uint uint32_t
#define FALSE 0
#define TRUE 1
#define CHAR_OFFSET 31
#define MIN_VALUE_BMAP2 87
#define DIST_UINT_MAX 4294967295ULL // max value of a 32-bit unsigned

typedef struct leaf {
	trie_uint count; // also tracks the cluster number
	trie_uint index;
	trie_uint edge_start; // start position of edges in the disk file
	dist_uint dist; // distance to original label, used for attentuation
} leaf;

typedef struct prefix {
	uint64_t bmap1 : 56; // 0, 32-86
	uint8_t count1 : 8;  // counts in this bitmap
	uint64_t bmap2 : 42; // 87-127
	uint8_t count2 : 8;   // counts in this bitmap
	void **child_nodes;  // 0 will be leaf, else will be prefix
} prefix;

prefix *initialize_trie(void);
trie_uint find_index_for_prefix_and_increment(char *s, prefix *trie, trie_uint* ctr, int should_increment);
void free_trie(prefix *trie);
#endif
