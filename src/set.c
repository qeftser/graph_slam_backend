
#include "set.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

set * construct_set() {
   /* allocate memory */
   set * ret = malloc(sizeof(set));
   ret->blocks = calloc(sizeof(set_block *),
                        SET_INITIAL_BLOCK_NUM);

   /* set initial values */
   ret->block_count = SET_INITIAL_BLOCK_NUM;

   return ret;
}

int insert_set(int element, set * s) {

   /* get the block this element is associated with */
   int block_num = element / SET_BLOCK_DIVISOR;

   if (block_num < s->block_count) {
insert_element:

      /* if our block has not been allocated,
       * allocate it                         */
      if (s->blocks[block_num] == NULL) {

         s->blocks[block_num] = calloc(sizeof(uint8_t),SET_BLOCK_SIZE);

         if (s->blocks[block_num] == NULL) {
            perror("[GRAPH SLAM BACKEND] failure allocating new set block");
            return -1; /* this is a bummer for the user */
         }
      }

      /* actually add the element */
      int block_pos = (element % SET_BLOCK_DIVISOR) / 8;
      int block_offset = element % 8;
      (*s->blocks[block_num])[block_pos] |= (1 << block_offset);

   }
   /* we have ran over our bounds for the array. We 
    * need to allocate more slots                  */
   else {
      do {

         /* double the number of slots */
         s->block_count *= 2;
         s->blocks = reallocarray(s->blocks,s->block_count,sizeof(set_block *));

         if (s->blocks == NULL) {
            perror("[GRAPH SLAM BACKEND] failure allocating more set slots");
            return -1; /* also a bummer for the user */
         }

         /* clear the new section */
         bzero((s->blocks + (s->block_count / 2)),sizeof(set_block *)*(s->block_count / 2));

      } while (block_num >= s->block_count); 
      goto insert_element;
   }

   return 0;
}

int member_set(int element, set * s) {

   /* get the block this element is associated with */
   int block_num = element / SET_BLOCK_DIVISOR;

   if (block_num < s->block_count && s->blocks[block_num] != NULL) {

      /* return the value in the element's position. It will
       * only be 1 if it has been added                     */
      int block_pos = (element % SET_BLOCK_DIVISOR) / 8;
      int block_offset = element % 8;
      return (*s->blocks[block_num])[block_pos] & (1 << block_offset);

   }

   /* Either our slot has not been allocated or our 
    * block has not been allocated. Either way, our
    * element is not in the set                     */
   return 0;
}


void remove_set(int element, set * s) {

   /* get the block this element is associated with */
   int block_num = element / SET_BLOCK_DIVISOR;

   if (block_num < s->block_count && s->blocks[block_num] != NULL) {

      /* mask out the element's position, removing it if it was set */
      int block_pos = (element % SET_BLOCK_DIVISOR) / 8;
      int block_offset = element % 8;
      (*s->blocks[block_num])[block_pos] &= (~(1 << block_offset));

   }

   /* Either our slot has not been allocated or our 
    * block has not been allocated. Either way, our
    * element is not in the set                     */
}

void destruct_set(set * s) {
   /* loop through all blocks, freeing them */
   for (int i = 0; i < s->block_count; ++i)
      free(s->blocks[i]); /* note that free(NULL) is allowed */

   /* free the array of blocks and the set */
   free(s->blocks);
   free(s);
}
