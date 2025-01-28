
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
   ret->block_num = 0;
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
      int block_pos = (element % SET_BLOCK_DIVISOR) / SET_BLOCK_SIZE;
      int block_offset = element % SET_BLOCK_SIZE;
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
