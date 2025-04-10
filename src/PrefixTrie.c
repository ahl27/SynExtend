#include "PrefixTrie.h"

leaf *alloc_leaf() {
	leaf *node = safe_malloc(sizeof(leaf));
	node->count = 0;
	node->index = 0;
	node->edge_start = 0;
	node->dist = 0;
	return node;
}

prefix *alloc_prefix() {
	prefix *p = safe_malloc(sizeof(prefix));
	p->bmap1 = 0;
	p->bmap2 = 0;
	p->count1 = 0;
	p->count2 = 0;
	p->child_nodes = NULL;
	return p;
}

prefix *initialize_trie(){
	return alloc_prefix();
}

void *insert_into_node_terminal(prefix *node){
	if(!(node->bmap1 & 1)){
		// create a leaf if this node wasn't previously terminal
		uint8_t total_children = node->count1 + node->count2;
		node->bmap1 |= 1;
		node->count1++;
		leaf *newleaf = alloc_leaf();
		void **new_ptr_holder = safe_malloc(sizeof(void*)*(total_children+1));
		new_ptr_holder[0] = newleaf;
		if(total_children){
			memcpy(&(new_ptr_holder[1]), node->child_nodes, sizeof(void*)*total_children);
			free(node->child_nodes);
		}
		node->child_nodes = new_ptr_holder;
	}
	return node->child_nodes[0];
}

void *insert_into_node_nonterminal(prefix *node, char s){
	uint8_t USE_HIGHER = s >= MIN_VALUE_BMAP2;
	// idx is either s-31 or s-86 depending on if we're using
	// the higher or lower bitfield
	uint8_t idx = s-(!USE_HIGHER*CHAR_OFFSET)-(USE_HIGHER*MIN_VALUE_BMAP2);
	uint8_t current_bit = 0;
	uint8_t insert_point = 0;
	uint64_t bitfield = node->bmap1;

	if(USE_HIGHER){
		insert_point = node->count1;
		bitfield = node->bmap2;
	}


	while(current_bit != idx){
		current_bit++;
		insert_point += bitfield & 1;
		bitfield >>= 1;
	}

	if(bitfield & 1){
		return node->child_nodes[insert_point];
	} else {
		uint8_t total_children = node->count1 + node->count2;
		// got to current bit BUT it doesn't exist

		// first create a new holder with an extra space
		void **new_ptr_holder = safe_malloc(sizeof(void*)*(total_children + 1));

		// then copy first <insert_point> pointers
		if(insert_point)
			memcpy(new_ptr_holder, node->child_nodes, sizeof(void*)*insert_point);

		// then add in the new pointer
		prefix *new_child = alloc_prefix();
		new_ptr_holder[insert_point] = new_child;

		// then copy the remaining pointers
		total_children -= insert_point;
		if(total_children)
			memcpy(&(new_ptr_holder[insert_point+1]),
							&(node->child_nodes[insert_point]),
							sizeof(void*)*total_children);

		// reassign the child pointer array
		if(node->child_nodes)
			free(node->child_nodes);
		node->child_nodes = new_ptr_holder;

		// trying to avoid weirdness from potential auto-cast to 32bit
		bitfield = 1;
		bitfield <<= idx;
		if(USE_HIGHER){
			node->count2++;
			node->bmap2 |= bitfield;
		} else {
			node->count1++;
			node->bmap1 |= bitfield;
		}
		return new_child;
	}
}


leaf *find_node_for_prefix(char *s, prefix *trie){
	prefix *tmp = trie;
	while(*s){
		if(*s < 31){
			free_trie(trie);
			error("Labels must contain ASCII values in the range 32-127 (received %u)", (uint8_t)(*s));
		}
		tmp = (prefix *)insert_into_node_nonterminal(tmp, *s++);
	}
	return (leaf *)insert_into_node_terminal(tmp);
}

trie_uint insert_into_trie(char *s, prefix *trie, trie_uint ctr, trie_uint to_add) {
	leaf *node = find_node_for_prefix(s, trie);
	if(!node->count)
		node->index = ctr++;
	node->count += to_add;
	return ctr;
}

trie_uint find_index_for_prefix(char *s, prefix *trie){
	leaf *l = find_node_for_prefix(s, trie);
	return l->index;
}

trie_uint find_index_for_prefix_and_increment(char *s, prefix *trie, trie_uint* ctr, int should_increment){
	// this function does three things:
	// 1. insert string into trie if it doesn't exist
	leaf *node = find_node_for_prefix(s, trie);
	if(!node->edge_start){
		node->index = (*ctr)++;
		node->edge_start = 1;
	}
	// 2. increment the in-degree of the node if incoming edge
	if(should_increment) node->count++;
	// 3. return the node's index
	return node->index;
}

void free_trie(prefix *trie) {
	if(!trie) return;

	uint8_t ctr = 0;
	uint8_t bits_remaining = trie->count1 + trie->count2;

	if(trie->bmap1 & 1)
		free((leaf*)(trie->child_nodes[ctr++]));

	while(ctr < bits_remaining)
		free_trie(trie->child_nodes[ctr++]);
	if(trie->child_nodes)
		free(trie->child_nodes);
	free(trie);
	return;
}
