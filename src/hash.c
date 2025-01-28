
#include "hash.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

hash * construct_table(int size) {

   /* allocate the table and fill standard values */
   hash * ret = malloc(sizeof(hash));
   ret->table_size = size;
   ret->entry = 0;
   ret->entry_ceiling = size;
   
   /* start by assuming there will only be size
    * entries, a load factor of 1              */
   ret->table = malloc(sizeof(int)*size);
   ret->entries = malloc(sizeof(he)*size);

   /* set all table values to HASH_SLOT_EMPTY */
   memset(ret->table,0xff,sizeof(int)*size);
   return ret;
}

void add_entry(long key, int val, hash * table) {

   /* perform FNV-1a hash on key */
   unsigned long hash = 0xcbf29ce484222325;
   int key2 = key;

   /* assume length of 8 bytes for long value */
   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;


   /* get table position */
   unsigned int pos = table->table[(hash % table->table_size)];
   /* check if slot is empty */
   if (pos == HASH_SLOT_EMPTY) {
      /* allocate a new value */
      pos = table->entry;

      /* register new element */
      table->table[(hash % table->table_size)] = pos;
      goto fill_slot;
   }
   /* slot is not empty, recurse through the hash
    * chain until we either find our key or hit an end */
   else {
      while (1) {

         /* is this our value? */
         if (table->entries[pos].key == key) {
            table->entries[pos].val = val;
            return; /* we are done */
         }
         /* is this the end of the chain? */
         else if (table->entries[pos].nextptr == HASH_CHAIN_END) {
            /* register the element */
            table->entries[pos].nextptr = table->entry;
            /* allocate the value */
            pos = table->entry;
            goto fill_slot;
         }

         /* move forward in the chain */
         pos = table->entries[pos].nextptr;
      }
   }

   /* fill the allocated slot with our values */
fill_slot:
   {
      table->entries[pos].key = key;
      table->entries[pos].val = val;
      table->entries[pos].nextptr = HASH_CHAIN_END;

      /* check if we have allocated the total number of entries
       * and allocate more entries if so                        */
      if (++(table->entry) == table->entry_ceiling) {
         table->entry_ceiling *= 2;
         table->entries = reallocarray(table->entries,table->entry_ceiling,sizeof(he));
         if (table->entries == NULL) {
            perror("[GRAPH SLAM BACKEND] hash table entry reallocation failed");
            return; /* crash is imminent */
         }
      }
   }
}

int get_entry(long key, hash * table) {

   /* perform FNV-1a hash on key */
   unsigned long hash = 0xcbf29ce484222325;
   int key2 = key;

   /* assume length of 8 bytes for long value */
   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;
   key2 >>= 8;

   hash ^= (key2 & 0xff);
   hash *= 0x00000100000001b3;

   /* get table position */
   unsigned int pos = table->table[(hash % table->table_size)];
   /* check if slot is empty */
   if (pos == HASH_SLOT_EMPTY) {

      /* the value does not exist, return failure */
      return -1;

   }
   /* slot is not empty, recurse through the hash
    * chain until we either find our key or hit an end */
   else {
      while (1) {

         /* is this our value? */
         if (table->entries[pos].key == key) {
            return table->entries[pos].val; /* we are done */
         }

         /* is this the end of the chain? */
         else if (table->entries[pos].nextptr == HASH_CHAIN_END) {

            /* the value does not exist, return failure */
            return -1;

         }

         /* move forward in the chain */
         pos = table->entries[pos].nextptr;
      }
   }
}

void clear_table(hash * table) {
   /* reset entry number */
   table->entry = 0;

   /* set all table slots to HASH_SLOT_EMPTY */
   memset(table->table,0xff,table->table_size*sizeof(int));
}

void destruct_table(hash * table) {
   /* free internal arrays */
   free(table->table);
   free(table->entries);
   /* free the table itself */
   free(table);
}
