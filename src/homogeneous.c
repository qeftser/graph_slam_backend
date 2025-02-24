
#include "homogeneous.h"
#include <math.h>
#include <stdio.h>

hv construct_homogeneous_vector(gsb_float x, gsb_float y) {
   hv ret = { x, y, 1 };
   return ret;
}

ht as_homogeneous_transformation(pv pos) {
   ht ret = { cos(pos.t), -sin(pos.t), pos.x,
              sin(pos.t),  cos(pos.t), pos.y,
              0,           0,          1     };
   return ret;
}

ht construct_homogeneous_translation(gsb_float trans_x, gsb_float trans_y) {
   ht ret = { 1, 0, trans_x,
              0, 1, trans_y,
              0, 0, 1       };
   return ret;
}

ht construct_homogeneous_rotation(gsb_float rad) {
   ht ret = { cos(rad), -sin(rad), 0,
              sin(rad),  cos(rad), 0,
              0,         0,        1 };
   return ret;
}

ht construct_homogeneous_transformation(gsb_float rad, gsb_float trans_x, gsb_float trans_y) {
   ht ret = { cos(rad), -sin(rad), trans_x,
              sin(rad),  cos(rad), trans_y,
              0,         0,        1       };
   return ret;
}

ht risky_merge_homogeneous_transformation(ht * t1, ht * t2) {
   ht ret = { t1->c1 * t2->c1 + t1->s1 * t2->s2,
              t1->c1 * t2->s1 + t1->s1 * t2->c2,
              t1->c1 * t2->tx + t1->s1 * t2->ty + t1->tx * t2->w,
              t1->s2 * t2->c1 + t1->c2 * t2->s2,
              t1->s2 * t2->s1 + t1->c2 * t2->c2,
              t1->s2 * t2->tx + t1->c2 * t2->ty + t1->ty * t2->w,
              0,
              0,
              t1->w * t2->w                                      };
   return ret;
}

ht merge_homogeneous_transformation(ht * t1, ht * t2) {
   ht ret = { t1->v[0] * t2->v[0] + t1->v[1] * t2->v[3] + t1->v[2] * t2->v[6],
              t1->v[0] * t2->v[1] + t1->v[1] * t2->v[4] + t1->v[2] * t2->v[7],
              t1->v[0] * t2->v[2] + t1->v[1] * t2->v[5] + t1->v[2] * t2->v[8],
              t1->v[3] * t2->v[0] + t1->v[4] * t2->v[3] + t1->v[5] * t2->v[6],
              t1->v[3] * t2->v[1] + t1->v[4] * t2->v[4] + t1->v[5] * t2->v[7],
              t1->v[3] * t2->v[2] + t1->v[4] * t2->v[5] + t1->v[5] * t2->v[8],
              t1->v[6] * t2->v[0] + t1->v[7] * t2->v[3] + t1->v[8] * t2->v[6],
              t1->v[6] * t2->v[1] + t1->v[7] * t2->v[4] + t1->v[8] * t2->v[7],
              t1->v[6] * t2->v[2] + t1->v[7] * t2->v[5] + t1->v[8] * t2->v[8] };
   return ret;
}

ht invert_homogeneous_transformation(ht * t) {
   ht ret = {  t->v[4] * t->v[8] - t->v[5] * t->v[7],
             -(t->v[1] * t->v[8] - t->v[2] * t->v[7]),
               t->v[1] * t->v[5] - t->v[2] * t->v[4],
             -(t->v[3] * t->v[8] - t->v[5] * t->v[6]),
               t->v[0] * t->v[8] - t->v[2] * t->v[6],
             -(t->v[0] * t->v[5] - t->v[2] * t->v[3]),
               t->v[3] * t->v[7] - t->v[4] * t->v[6],
             -(t->v[0] * t->v[7] - t->v[1] * t->v[6]),
               t->v[0] * t->v[4] - t->v[1] * t->v[3]  };
   gsb_float det = t->v[0] * ret.v[0] + t->v[3] * ret.v[3] + t->v[6] * ret.v[6];
   ret.v[0] /= det; ret.v[1] /= det; ret.v[2] /= det;
   ret.v[3] /= det; ret.v[4] /= det; ret.v[5] /= det;
   ret.v[6] /= det; ret.v[7] /= det; ret.v[8] /= det;
   return ret;
}

hv apply_homogeneous_transformation(ht * t, hv * p) {
   hv ret = { t->c1 * p->x + t->s1 * p->y + t->tx * p->w,
              t->s2 * p->x + t->c2 * p->y + t->ty * p->w,
              t->w  * p->w                               };
   return ret;
}

lv destruct_homogeneous_vector(hv vector) {
   lv ret = { vector.x / vector.w, vector.y / vector.w };
   return ret;
}

pv destruct_homogeneous_transformation(ht * t) {
   pv ret = { t->tx / t->w, t->ty / t->w, acos(t->c1 / t->w) };
   return ret;
}

void normalize_homogeneous_transformation(ht * t) {
   t->v[0] /= t->v[8]; t->v[1] /= t->v[8]; t->v[2] /= t->v[8];
   t->v[3] /= t->v[8]; t->v[4] /= t->v[8]; t->v[5] /= t->v[8];
   t->v[6] /= t->v[8]; t->v[7] /= t->v[8]; t->v[8] = 1;
}

void normalize_homogeneous_vector(hv * v) {
   v->x /= v->w; v->y /= v->w; v->w = 1;
}

void print_homogeneous_transformation(ht * t) {
   printf("transformation values:\n");
   printf("[ %02.3f %02.3f %02.3f ]\n",t->v[0],t->v[1],t->v[2]);
   printf("[ %02.3f %02.3f %02.3f ]\n",t->v[3],t->v[4],t->v[5]);
   printf("[ %02.3f %02.3f %02.3f ]\n",t->v[6],t->v[7],t->v[8]);
   printf("equivilant euclidean values:\n");
   pv euclidean = destruct_homogeneous_transformation(t);
   printf("x: %02.3f y: %02.3f 0: %02.3f\n",
         euclidean.x,euclidean.y,euclidean.t);
}

void print_homogeneous_vector(hv * v) {
   printf("[ %02.3f %02.3f %02.3f ]\n",v->x,v->y,v->w);
}
