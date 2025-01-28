
#ifndef DEFAULT

#define DEFAULT

/* holds the x, y and theta values
 * of a position. Can also be seen
 * as a 1x3 vector                 */
typedef struct position_vector {
   float x;
   float y;
   float t;
} pv;

/* cartesian point */
typedef struct location_vector {
   float x;
   float y;
} lv;

#endif
