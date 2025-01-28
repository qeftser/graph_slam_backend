
#ifndef HASH

/* Simple long->int hash table that keeps track of
 * entries in the H matrix for the construction phase.
 *
 * This hash table is intended to be used with the pattern
 * load -> clear -> load
 * As a result of this it does not do anything to prevent
 * memory fragmentation. Therefore, it should not be
 * used when the pattern calls for lots of removals.
 */
#define HASH

/* Internal macro. Indicates that the given
 * hash slot is empty                      */
#define HASH_SLOT_EMPTY 0xffffffff
/* Internal macro. Indicates that the given
 * hash entry is the end of a chain        */
#define HASH_CHAIN_END  -1

/* An entry in the table. The key is hashed
 * so as to avoid performance drops with
 * linear keys. nextptr is actually an offset
 * in the entries array of the hash_table
 * struct. This is done to keep the entries
 * small and to allow bulk allocation of nodes
 * to reduce calls to malloc                  */
typedef struct hash_entry {
   long key;
   int  val;
   int  nextptr;
} he;

/* The hash table. The table array
 * holds a list of offsets into the 
 * entries array, which holds the
 * hash entries. There are table_size
 * table slots, and entry_ceiling
 * entries slots. entry keeps track
 * of the next entry to allocate. It is
 * basically like a stack pointer       */
typedef struct hash_table {
   int table_size;
   int entry, entry_ceiling;
   int * table;
   he * entries;
} hash;

/* construct and return a new hash table of
 * the give size                           */
hash * construct_table(int size);
/* insert a new entry into the hash table. If
 * the given key already exists, replace the
 * existing value                            */
void add_entry(long key, int val, hash * table);
/* return the value stored with the given key.
 * return -1 if the key does not exist.       */
int  get_entry(long key, hash * table);
/* remove all entries from the hash table */
void clear_table(hash * table);
/* destroy and free all memory associated with
 * the hash table. Note that this also frees
 * the memory held by the passed pointer      */
void destruct_table(hash * table);




#endif
