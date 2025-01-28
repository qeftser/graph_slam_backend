
#ifndef SET

/* A set. This is for simple membership checking
 * and is used in the graph slam algorithm to
 * check whether or not a node is fixed in the
 * optimization. The implimentation is a bit array
 * that is broken into chuncks and referenced by
 * another array. The first array is grown as more
 * chunks are needed and the second array is allocated
 * when an element is accessed.                       */

/* Note that this set is only really effective for elements
 * that are expected to be accessed in a linear pattern, i.e.
 * no jumping from element 1 to 200,000. This is because
 * accesses are done as array accesses, so there is no real
 * good way to efficently add values that are very far away
 * from each other                                           */

#define SET
#include <stdint.h>

/* each block will have 64 bytes in it, which
 * translates to 512 elements in a block.    */
#define SET_BLOCK_SIZE 64
#define SET_BLOCK_DIVISOR (SET_BLOCK_SIZE*sizeof(uint8_t)*8)

/* start with 32 slots in our array. This is
 * 16384 elements                            */
#define SET_INITIAL_BLOCK_NUM 32

/* a block in our set */
typedef uint8_t set_block[SET_BLOCK_SIZE];

/* The set structure. An expanding array
 * of set blocks                         */
typedef struct bit_array_set {
   int block_num, block_count;
   set_block ** blocks;
} set;

/* allocate a new set */
set * construct_set(void);

/* add an element to our set */
int insert_set(int element, set * s);

/* remove an element from our set */
int remove_set(int element, set * s);

/* free all memory associated with our set.
 * Note that this also frees the reference
 * s                                      */
void destruct_set(set * s);




#endif
