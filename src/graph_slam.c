
#include "graph_slam.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* helper method for debugging */
void print_m33(m33 * m) {
   printf("[ %f %f %f ]\n",(*m)[0],(*m)[1],(*m)[2]);
   printf("[ %f %f %f ]\n",(*m)[3],(*m)[4],(*m)[5]);
   printf("[ %f %f %f ]\n",(*m)[6],(*m)[7],(*m)[8]);
}

void compute_jacobians_and_error(pge * edge, pg * graph, m33 * J_i, m33 * J_j, pv * e_ij) {
   pv x_i = graph->node[edge->xi].pos;
   pv x_j = graph->node[edge->xj].pos;
   double cos_iij = cos(edge->observation.t + x_i.t);
   double sin_iij = sin(edge->observation.t + x_i.t);
   double cos_ij = cos(edge->observation.t);
   double sin_ij = sin(edge->observation.t);
   double cos_i = cos(x_i.t);
   double sin_i = sin(x_i.t);

   (*e_ij).x = (x_j.x - x_i.x)*cos_iij + (x_j.y - x_i.y)*sin_iij 
             - edge->observation.y*sin_ij - edge->observation.x*cos_ij;
   (*e_ij).y = (x_i.x - x_j.x)*sin_iij + (x_j.y - x_i.y)*cos_iij
             + edge->observation.x*sin_ij - edge->observation.y*cos_ij;
   (*e_ij).t = x_j.t - x_i.t - edge->observation.t;
   /* normalize angle. We get bad results otherwise */
   normalize_angle(e_ij->t);

   /* de_ij / dx_i */
   (*J_i)[0] = -cos_iij;
   (*J_i)[3] =  sin_iij;
   (*J_i)[6] =  0;

   /* de_ij / dy_i */
   (*J_i)[1] = -sin_iij;
   (*J_i)[4] = -cos_iij;
   (*J_i)[7] =  0;

   /* de_ij / d0_i */
   (*J_i)[2] = (x_i.x - x_j.x)*sin_i + (x_j.y - x_i.y)*cos_i;
   (*J_i)[5] = (x_i.x - x_j.x)*cos_i + (x_i.y - x_j.y)*sin_i;
   (*J_i)[8] = -1;

   /* de_ij / dx_j */
   (*J_j)[0] =  cos_iij;
   (*J_j)[3] = -sin_iij;
   (*J_j)[6] =  0;

   /* de_ij / dy_j */
   (*J_j)[1] =  sin_iij;
   (*J_j)[4] =  cos_iij;
   (*J_j)[7] =  0;

   /* de_ij / d0_j */
   (*J_j)[2] =  0;
   (*J_j)[5] =  0;
   (*J_j)[8] =  1;

}

void compute_constraints_and_coefficients
(m33 * J_i, m33 * J_j, pv * e_ij, sm33 * i, m33 * H_ii, m33 * H_ji, m33 * H_jj, pv * b_i, pv * b_j) {

   m33 temp;

   /* H_ii = (J_i)' i J_i */
   temp[0] = (*J_i)[0]*i->a00 + (*J_i)[3]*i->a10 + (*J_i)[6]*i->a20;
   temp[1] = (*J_i)[0]*i->a10 + (*J_i)[3]*i->a11 + (*J_i)[6]*i->a21;
   temp[2] = (*J_i)[0]*i->a20 + (*J_i)[3]*i->a21 + (*J_i)[6]*i->a22;
   temp[3] = (*J_i)[1]*i->a00 + (*J_i)[4]*i->a10 + (*J_i)[7]*i->a20;
   temp[4] = (*J_i)[1]*i->a10 + (*J_i)[4]*i->a11 + (*J_i)[7]*i->a21;
   temp[5] = (*J_i)[1]*i->a20 + (*J_i)[4]*i->a21 + (*J_i)[7]*i->a22;
   temp[6] = (*J_i)[2]*i->a00 + (*J_i)[5]*i->a10 + (*J_i)[8]*i->a20;
   temp[7] = (*J_i)[2]*i->a10 + (*J_i)[5]*i->a11 + (*J_i)[8]*i->a21;
   temp[8] = (*J_i)[2]*i->a20 + (*J_i)[5]*i->a21 + (*J_i)[8]*i->a22;
   (*H_ii)[0] = temp[0]*(*J_i)[0] + temp[1]*(*J_i)[3] + temp[2]*(*J_i)[6];
   (*H_ii)[1] = temp[0]*(*J_i)[1] + temp[1]*(*J_i)[4] + temp[2]*(*J_i)[7];
   (*H_ii)[2] = temp[0]*(*J_i)[2] + temp[1]*(*J_i)[5] + temp[2]*(*J_i)[8];
   (*H_ii)[3] = temp[3]*(*J_i)[0] + temp[4]*(*J_i)[3] + temp[5]*(*J_i)[6];
   (*H_ii)[4] = temp[3]*(*J_i)[1] + temp[4]*(*J_i)[4] + temp[5]*(*J_i)[7];
   (*H_ii)[5] = temp[3]*(*J_i)[2] + temp[4]*(*J_i)[5] + temp[5]*(*J_i)[8];
   (*H_ii)[6] = temp[6]*(*J_i)[0] + temp[7]*(*J_i)[3] + temp[8]*(*J_i)[6];
   (*H_ii)[7] = temp[6]*(*J_i)[1] + temp[7]*(*J_i)[4] + temp[8]*(*J_i)[7];
   (*H_ii)[8] = temp[6]*(*J_i)[2] + temp[7]*(*J_i)[5] + temp[8]*(*J_i)[8];


   /* b_i = (J_i)' i e_ij */
   b_i->x = temp[0]*e_ij->x + temp[1]*e_ij->y + temp[2]*e_ij->t;
   b_i->y = temp[3]*e_ij->x + temp[4]*e_ij->y + temp[5]*e_ij->t;
   b_i->t = temp[6]*e_ij->x + temp[7]*e_ij->y + temp[8]*e_ij->t;

   /* H_jj = (J_j)' i J_j */
   temp[0] = (*J_j)[0]*i->a00 + (*J_j)[3]*i->a10 + (*J_j)[6]*i->a20;
   temp[1] = (*J_j)[0]*i->a10 + (*J_j)[3]*i->a11 + (*J_j)[6]*i->a21;
   temp[2] = (*J_j)[0]*i->a20 + (*J_j)[3]*i->a21 + (*J_j)[6]*i->a22;
   temp[3] = (*J_j)[1]*i->a00 + (*J_j)[4]*i->a10 + (*J_j)[7]*i->a20;
   temp[4] = (*J_j)[1]*i->a10 + (*J_j)[4]*i->a11 + (*J_j)[7]*i->a21;
   temp[5] = (*J_j)[1]*i->a20 + (*J_j)[4]*i->a21 + (*J_j)[7]*i->a22;
   temp[6] = (*J_j)[2]*i->a00 + (*J_j)[5]*i->a10 + (*J_j)[8]*i->a20;
   temp[7] = (*J_j)[2]*i->a10 + (*J_j)[5]*i->a11 + (*J_j)[8]*i->a21;
   temp[8] = (*J_j)[2]*i->a20 + (*J_j)[5]*i->a21 + (*J_j)[8]*i->a22;
   (*H_jj)[0] = temp[0]*(*J_j)[0] + temp[1]*(*J_j)[3] + temp[2]*(*J_j)[6];
   (*H_jj)[1] = temp[0]*(*J_j)[1] + temp[1]*(*J_j)[4] + temp[2]*(*J_j)[7];
   (*H_jj)[2] = temp[0]*(*J_j)[2] + temp[1]*(*J_j)[5] + temp[2]*(*J_j)[8];
   (*H_jj)[3] = temp[3]*(*J_j)[0] + temp[4]*(*J_j)[3] + temp[5]*(*J_j)[6];
   (*H_jj)[4] = temp[3]*(*J_j)[1] + temp[4]*(*J_j)[4] + temp[5]*(*J_j)[7];
   (*H_jj)[5] = temp[3]*(*J_j)[2] + temp[4]*(*J_j)[5] + temp[5]*(*J_j)[8];
   (*H_jj)[6] = temp[6]*(*J_j)[0] + temp[7]*(*J_j)[3] + temp[8]*(*J_j)[6];
   (*H_jj)[7] = temp[6]*(*J_j)[1] + temp[7]*(*J_j)[4] + temp[8]*(*J_j)[7];
   (*H_jj)[8] = temp[6]*(*J_j)[2] + temp[7]*(*J_j)[5] + temp[8]*(*J_j)[8];

   /* H_ji = (J_j)' i J_i */
   (*H_ji)[0] = temp[0]*(*J_i)[0] + temp[1]*(*J_i)[3] + temp[2]*(*J_i)[6];
   (*H_ji)[1] = temp[0]*(*J_i)[1] + temp[1]*(*J_i)[4] + temp[2]*(*J_i)[7];
   (*H_ji)[2] = temp[0]*(*J_i)[2] + temp[1]*(*J_i)[5] + temp[2]*(*J_i)[8];
   (*H_ji)[3] = temp[3]*(*J_i)[0] + temp[4]*(*J_i)[3] + temp[5]*(*J_i)[6];
   (*H_ji)[4] = temp[3]*(*J_i)[1] + temp[4]*(*J_i)[4] + temp[5]*(*J_i)[7];
   (*H_ji)[5] = temp[3]*(*J_i)[2] + temp[4]*(*J_i)[5] + temp[5]*(*J_i)[8];
   (*H_ji)[6] = temp[6]*(*J_i)[0] + temp[7]*(*J_i)[3] + temp[8]*(*J_i)[6];
   (*H_ji)[7] = temp[6]*(*J_i)[1] + temp[7]*(*J_i)[4] + temp[8]*(*J_i)[7];
   (*H_ji)[8] = temp[6]*(*J_i)[2] + temp[7]*(*J_i)[5] + temp[8]*(*J_i)[8];

   /* b_j = (J_j)' i e_ij */
   b_j->x = temp[0]*e_ij->x + temp[1]*e_ij->y + temp[2]*e_ij->t;
   b_j->y = temp[3]*e_ij->x + temp[4]*e_ij->y + temp[5]*e_ij->t;
   b_j->t = temp[6]*e_ij->x + temp[7]*e_ij->y + temp[8]*e_ij->t;

}

int insert_information_matrix_values(int xi, int xj, m33 * H_ii, m33 * H_ji, m33 * H_jj, mt * values, hash * table, int value) {
   int v_pos, h_pos_i = xi * 3, h_pos_j = xj * 3;
   /* insert the 3x3 matrix H_ii into the list of triplets. Note
    * that we need to check for and mark existance of these
    * values as we may have to add to them later on and cannot
    * explicitly allocate the matrix.                           */
   if ( (v_pos = get_entry(long_from_ints(xi,xi),table)) == -1) {
      /* there is no curent entry for this value, make one */
      v_pos = value;
      value += 6; /* only 6 values for H_ii, it is on the
                     diagonal, so the upper half is symmetric */
      add_entry(long_from_ints(xi,xi),v_pos,table);

      /* set the column and row indices for these values
       * in the final matrix H. These values are laid out
       * the same way as the symmetric3x3 struct. We proceed
       * down a row and then drop rows when we hit the diagonal
       * value                                                 */
      values[v_pos  ].col = h_pos_i  ; values[v_pos  ].row = h_pos_i  ;
      values[v_pos+1].col = h_pos_i  ; values[v_pos+1].row = h_pos_i+1;
      values[v_pos+2].col = h_pos_i+1; values[v_pos+2].row = h_pos_i+1;
      values[v_pos+3].col = h_pos_i  ; values[v_pos+3].row = h_pos_i+2;
      values[v_pos+4].col = h_pos_i+1; values[v_pos+4].row = h_pos_i+2;
      values[v_pos+5].col = h_pos_i+2; values[v_pos+5].row = h_pos_i+2;
   }
   /* now we can insert the values. We will go by 
    * rows, so H_xx[1][2] = H_xx[4]*/
   values[v_pos  ].val += (*H_ii)[0];
   values[v_pos+1].val += (*H_ii)[3];
   values[v_pos+2].val += (*H_ii)[4];
   values[v_pos+3].val += (*H_ii)[6];
   values[v_pos+4].val += (*H_ii)[7];
   values[v_pos+5].val += (*H_ii)[8];

   /* repeat our actions for the matrix H_jj */
   if ( (v_pos = get_entry(long_from_ints(xj,xj),table)) == -1) {
      v_pos = value;
      value += 6; 
      add_entry(long_from_ints(xj,xj),v_pos,table);

      values[v_pos  ].col = h_pos_j  ; values[v_pos  ].row = h_pos_j  ;
      values[v_pos+1].col = h_pos_j  ; values[v_pos+1].row = h_pos_j+1;
      values[v_pos+2].col = h_pos_j+1; values[v_pos+2].row = h_pos_j+1;
      values[v_pos+3].col = h_pos_j  ; values[v_pos+3].row = h_pos_j+2;
      values[v_pos+4].col = h_pos_j+1; values[v_pos+4].row = h_pos_j+2;
      values[v_pos+5].col = h_pos_j+2; values[v_pos+5].row = h_pos_j+2;
   }
   values[v_pos  ].val += (*H_jj)[0];
   values[v_pos+1].val += (*H_jj)[3];
   values[v_pos+2].val += (*H_jj)[4];
   values[v_pos+3].val += (*H_jj)[6];
   values[v_pos+4].val += (*H_jj)[7];
   values[v_pos+5].val += (*H_jj)[8];

   /* now we proceed to the matrix H_ji. This
    * one is a little different, as it not symmetric
    * and should sit below the diagonal. Because of
    * this, we need to swap xj and xi if xj > xi. This
    * will work because the matrix we are passed as H_ji
    * is assumed to be the lower diagonal one. If this is
    * not the case, this assumption is incorect and the 
    * information matrix will be ill-formed              */
   if (xi > xj) {
      int temp = xj;
      xj = xi;
      xi = temp;
      temp = h_pos_i;
      h_pos_i = h_pos_j;
      h_pos_j = temp;

   }
   /* We can proceed in an identical fashion as the previous blocks,
    * except that we insert the entire matrix instead of just the diagonal
    */
   if ( (v_pos = get_entry(long_from_ints(xi,xj),table)) == -1) {
      v_pos = value;
      value += 9; /* we do 9 this time because we are adding a 
                     full 3x3 matrix                          */
      add_entry(long_from_ints(xi,xj),v_pos,table);

      values[v_pos  ].col = h_pos_i  ; values[v_pos  ].row = h_pos_j  ;
      values[v_pos+1].col = h_pos_i+1; values[v_pos+1].row = h_pos_j  ;
      values[v_pos+2].col = h_pos_i+2; values[v_pos+2].row = h_pos_j  ;
      values[v_pos+3].col = h_pos_i  ; values[v_pos+3].row = h_pos_j+1;
      values[v_pos+4].col = h_pos_i+1; values[v_pos+4].row = h_pos_j+1;
      values[v_pos+5].col = h_pos_i+2; values[v_pos+5].row = h_pos_j+1;
      values[v_pos+6].col = h_pos_i  ; values[v_pos+6].row = h_pos_j+2;
      values[v_pos+7].col = h_pos_i+1; values[v_pos+7].row = h_pos_j+2;
      values[v_pos+8].col = h_pos_i+2; values[v_pos+8].row = h_pos_j+2;
   }
   values[v_pos  ].val += (*H_ji)[0];
   values[v_pos+1].val += (*H_ji)[1];
   values[v_pos+2].val += (*H_ji)[2];
   values[v_pos+3].val += (*H_ji)[3];
   values[v_pos+4].val += (*H_ji)[4];
   values[v_pos+5].val += (*H_ji)[5];
   values[v_pos+6].val += (*H_ji)[6];
   values[v_pos+7].val += (*H_ji)[7];
   values[v_pos+8].val += (*H_ji)[8];

   /* we need to know to what position we
    * have advanced in our values vector */
   return value;
}

pgn_handle add_node(pv pos, pg * graph) {
   pgn_handle new = graph->node_count;

   /* check if we have hit the total allocated nodes.
    * If so, allocate more                           */
   if (++(graph->node_count) == graph->node_ceiling) {
      graph->node_ceiling *= 2;
      graph->node = reallocarray(graph->node,graph->node_ceiling,sizeof(pgn));
      if (graph->node == NULL) {
         perror("[GRAPH SLAM BACKEND] node array reallocation failed");
         return -1; /* this will caouse the system to crash shortly */
      }
   }
   /* load provided values into node */
   graph->node[new].pos = pos;
   return new;
}

int add_edge(pgn_handle xi, pgn_handle xj, pv * observation, sm33 * information, pg * graph) {
   int new = graph->edge_count;

   /* check if we have hit the total allocated edges.
    * If so, allocate more                           */
   if (++(graph->edge_count) == graph->edge_ceiling) {
      graph->edge_ceiling *= 2;
      graph->edge = reallocarray(graph->edge,graph->edge_ceiling,sizeof(pge));
      if (graph->edge == NULL) {
         perror("[GRAPH SLAM BACKEND] edge array reallocation failed");
         return -1; /* this will caouse the system to crash shortly */
      }
   }

   /* load provided values into edge */
   graph->edge[new].xi = xi;
   graph->edge[new].xj = xj;
   memcpy(&graph->edge[new].observation,observation,sizeof(pv));
   memcpy(&graph->edge[new].information,information,sizeof(sm33));
   return new;
}

pg * construct_pose_graph() {
   /* allocate our memory. We will start with 256
    * node and edge slots allocated              */
   pg * ret = malloc(sizeof(pg));
   ret->node = malloc(sizeof(pgn)*256);
   ret->edge = malloc(sizeof(pge)*256);

   /* set initial status values */
   ret->node_count = 0;
   ret->node_ceiling = 256;
   ret->edge_count = 0;
   ret->edge_ceiling = 256;

   return ret;
}

void free_pose_graph(pg * graph) {
   /* deallocate our memory and free
    * graph                         */
   free(graph->node);
   free(graph->edge);
   free(graph);
}

int optimize(pg * graph) {

   /* setup structures needed for computation */
   hash  * table = construct_table(graph->edge_count);
   float * b = malloc(sizeof(float)* 3 *graph->node_count);
   mt * values = malloc(sizeof(mt)* 21 *graph->edge_count);

   /* our matrices will be allocated when we determine the
    * number of non-zero elements in them                 */
   pcmat * H = NULL;
   pcmat * L = NULL;

   /* keep track of the total movement of the system */
   double total_change;
   double largest_change;

   /* perform optimization loop */
   while (1) {

      bzero(b,sizeof(float)* 3 *graph->node_count);
      bzero(values,sizeof(mt)* 21 *graph->edge_count);
      clear_table(table);
      total_change = 0.0;
      largest_change = 0.0;

      int value = 0, b_pos, xi, xj;
      m33 J_i, J_j, H_ii, H_ji, H_jj;
      pv e_ij, b_i, b_j;
      for (int i = 0; i < graph->edge_count; ++i) {

         xi = graph->edge[i].xi;
         xj = graph->edge[i].xj;

         /* perform computation to get our jacobian blocks H_ii, H_ji, H_jj
          * and our error values b_i and b_j. These are the parameters to the
          * equations we will use the choloskey decomposition to solve       */
         compute_jacobians_and_error(&graph->edge[i],graph,&J_i,&J_j,&e_ij);

         /* if xi > xj we need to swap our jacobians when we compute our
          * constraints. This will result in the production of the lower
          * diagonal H_ji as opposed to H_ij, keeping our values list to
          * only the lower diagonal ones                                 */
         if (xi > xj)
            compute_constraints_and_coefficients(&J_j,&J_i,&e_ij,&graph->edge[i].information,
                                                 &H_jj,&H_ji,&H_ii,&b_j,&b_i);
         else
            compute_constraints_and_coefficients(&J_i,&J_j,&e_ij,&graph->edge[i].information,
                                                 &H_ii,&H_ji,&H_jj,&b_i,&b_j);

         /* insert error at position i in the vector b 
          * note that each error contains three values,
          * x, y and 0. This is the reason for the
          * multiplication by three. Additionally, 
          * because in our final system we are solving
          * for -b, we will subtract these values instead
          * of adding them                               */
         b_pos = xi * 3;
         b[b_pos  ] -= b_i.x;
         b[b_pos+1] -= b_i.y;
         b[b_pos+2] -= b_i.t;

         /* insert error at position j in the vector b */
         b_pos = xj * 3;
         b[b_pos  ] -= b_j.x;
         b[b_pos+1] -= b_j.y;
         b[b_pos+2] -= b_j.t;

         /* insert matrices H_ii,H_ij,H_jj into the information matrix H */
         value = insert_information_matrix_values(xi,xj,&H_ii,&H_ji,&H_jj,values,table,value);
      }

      /* fix the first entry to make our system well defined. Because
       * our insertion function inserts the diagonal matricies first, 
       * we know that the first entry in values will be diagonal. Add
       * the identity matrix to this entry                           */
      values[0].val += 1;
      values[2].val += 1;
      values[5].val += 1;

      /* now we have enough information to allocate our matricies if
       * they have not yet been allocated. Note that value now holds
       * the total non-zero elements present in H                   */
      if (H == NULL) {
         H = construct_packed_column_matrix(values,graph->node_count*3,value);
         /* we can reuse this matrix for each remaining iteration */
         L = compute_symbolic_factorization(H); 
      }
      else
         load_packed_column_matrix(values,H);

      /* perform the factorization and solve the matrix */
      perform_numerical_factorization(H,L);
      solve_system(L,b);

      /* collect total and largest change */
      for (int i = 0; i < graph->node_count*3; ++i) {
         total_change += fabs(b[i]);
         if (fabs(b[i]) > largest_change)
            largest_change = fabs(b[i]);
      }

      /* check if our decomposition has diverged. This is done to avoid letting everything
       * go to nan. Unfortunatly we cannot tell if this will happen before it does, so the
       * user may be left with poses that are partially optimized and make no sense. There
       * is not much that can be done about that, the hope is that they will go fix their
       * front end if they see this message enough. Regardless, the guarentee this section
       * provides is that a call to optimize will never result in the poses becoming undefined */
      if (isnan(total_change) || isinf(total_change) || total_change != total_change) {
         printf("[GRAPH SLAM BACKEND] optimization failure, pose graph is likely not positive semi-definite\n");
         printf("[GRAPH SLAM BACKEND] exiting prematurely, pose graph values may be trash\n");
         return -1;
      }

      /* update our poses given the solution to the 
       * problem we constructed                    */
      for (int i = 0; i < graph->node_count; ++i) {
         graph->node[i].pos.x += b[i*3];
         graph->node[i].pos.y += b[(i*3)+1];
         graph->node[i].pos.t += b[(i*3)+2];
         normalize_angle(graph->node[i].pos.t);
      }

      /* temporary values. Print all updates and hold for 
       * user to press enter                             */
      printf("=========================================\n");
      for (int i = 0; i < graph->node_count; ++i) {
         printf("[ %.2f %.2f %.2f ] [ %.2f %.2f %.2f ]\n",graph->node[i].pos.x,graph->node[i].pos.y,
               graph->node[i].pos.t,b[i*3],b[(i*3)+1],b[(i*3)+2]);
      }
      printf("=========================================\n");
      printf("TOTAL CHANGE: %e\n",total_change);
      printf("LARGE CHANGE: %e\n",largest_change);
      printf("=========================================\n");
      getchar();
   }

   /* cleanup all the memory we allocated */
   destruct_table(table);
   free(b);
   free(values);
   free_packed_column_matrix(H);
   free_packed_column_matrix(L);

   return 0;

}
