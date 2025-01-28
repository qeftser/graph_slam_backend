
#ifndef HOMOGENEOUS

/* a simple set of functions for computing
 * two dimensional (3 dof) homogeneous
 * transformations and coordinates. This is
 * by no means a full fledged set of values
 * but does pretty good for what it is      */

/* references:
 * https://en.wikipedia.org/wiki/Homogeneous_coordinates
 * https://www.youtube.com/watch?v=MQdm0Z_gNcw&list=PLgnQpQtFTOGRYjqjdZxTEQPZuFHQa7O7Y&index=30
 */

#define HOMOGENEOUS
#include "default.h"

/* a vector in homogeneous representation.*/
typedef struct homogeneous_vector {
   float x;
   float y;
   float w;
} hv;

/* a matrix in homogeneous representation.
 * Note that this is just a 3x3 matrix. Note
 * also that there are two ways of accessing
 * the same values                           */
typedef struct homogeneous_transform {
   union {
      struct {
         float c1, s1, tx;
         float s2, c2, ty;
         float z1, z2, w;
      };
      float v[9];
   };
} ht;

/* return a homogeneous vector given x and y */
hv construct_homogeneous_vector(float x, float y);
/* return a homogeneous transformation given a 
 * 2d position vector                          */
ht as_homogeneous_transformation(pv pos);
/* return a homogeneous transformation with a rotation
 * of zero and x and y translations as given           */
ht construct_homogeneous_translation(float trans_x, float trans_y);
/* return a homogeneous transformation with a translation
 * of zero and a rotation of the given radians           */
ht construct_homogeneous_rotation(float rad);
/* produce a full homogeneous transformation with x and y tranlations
 * defined, as well the rotation in radians                          */
ht construct_homogeneous_transformation(float rad, float trans_x, float trans_y);
/* perform 3x3 matrix multiplication to combine two homogeneous
 * transformations, producing a new one                        */
ht merge_homogeneous_transformation(ht * t1, ht * t2);
/* perform 3x3 matrix inversion, returning the inverse of the
 * given homogeneous transformation                          */
ht invert_homogeneous_transformation(ht * t);
/* perform the given homogeneous transformation on the homogeneous
 * vector p. Equvilant to 3x3 * 3x1 matrix vector multiplication.
 * Return the resulting vector                                     */
hv apply_homogeneous_transformation(ht * t, hv * p);
/* recover the normal cartesian point given
 * a homogeneous vector. Returns [ x/w y/w ] */
lv destruct_homogeneous_vector(hv vector);
/* recover the normal 2d position given a homogeneous
 * transformation matrix. Performs [ ht[2]/w ht[5]/w acos(ht[0]) ] */
pv destruct_homogeneous_transformation(ht * t);
/* divides all values in the transformation by w */
void normalize_homogeneous_transformation(ht * t);
/* divides all values in the vector by w */
void normalize_homogeneous_vector(hv * v);
/* helper method to print a homogeneous transformation
 * Also prints the equivilant position vector, which is handy */
void print_homogeneous_transformation(ht * t);
/* helper method to print a homogeneous vector */
void print_homogeneous_vector(hv * v);


#endif
