
#ifndef DEFAULT

#define DEFAULT

typedef double gsb_float;

/* holds the x, y and theta values
 * of a position. Can also be seen
 * as a 1x3 vector                 */
typedef struct position_vector {
   gsb_float x;
   gsb_float y;
   gsb_float t;
} pv;

/* cartesian point */
typedef struct location_vector {
   gsb_float x;
   gsb_float y;
} lv;

#endif
