
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

/* sucessor function of the above two functions */
int compute_constraints_and_fill_information_matrix(int xi, int xj, m33 * J_i, m33 * J_j, pv * e_ij, sm33 * i,
                                                    float * b, mt * values, hash * table, set * fixed,
                                                    int value) {
   /* intermediate values */
   m33 temp, H_xx;
   sm33 H_ii, H_jj;
   pv b_i, b_j;

   /* values for keeping track of our position in the values 
    * array.                                                */
   int v_pos, h_pos_i = xi * 3, h_pos_j = xj * 3;

   /* we are doing some fancy bitwise stuff here but it is normally
    * the best way to do these things. The value 1 will indicate 
    * that xi is a fixed node and the value 2 will indicate that xj
    * is fixed. The value 3 will also occur, and this is when both
    * values are fixed. Here we simply exit because no values will
    * be added to our information matrix                            */
   int usage = (member_set(xi,fixed) | (member_set(xj,fixed) << 1));

   /* perform operations based on which values are fixed and which aren't */
   switch(usage) {
      /* no values are fixed. Compute all of H_ii, H_ji, H_jj and insert
       * it into the values array.                                      */
      case 0:
         /* Because we only store the lower diagonal portion of the matrix H,
          * we need to be careful about which block we compute. Our goal is
          * H_ji, but if xi is greater than xj we would end up computing
          * H_ij instead, which is in the upper diagonal area of H. To avoid
          * this, we will swap some of our operations and actually compute
          * H_ji for this case, hence the branch.                           */
         if (xi > xj) {

            /* H_ii = (J_i)' i J_i */
            m33_A_transpose_mult_B_symmetric(temp,(*J_i),(*i));
            sm33_A_mult_B(H_ii,temp,(*J_i));

            /* b_i = (J_i)' i e_ij */
            pv_A_transpose_mult_b(b_i,temp,(*e_ij));

            /* H_ij = (J_i)' i J_j */
            m33_A_mult_B(H_xx,temp,(*J_j));

            /* H_jj = (J_j)' i J_j */
            m33_A_transpose_mult_B_symmetric(temp,(*J_j),(*i));
            sm33_A_mult_B(H_jj,temp,(*J_j));

            /* b_j = (J_j)' i e_ij */
            pv_A_transpose_mult_b(b_j,temp,(*e_ij));

            /* Insert the 3x3 matrix H_ji into the values array, allocating
             * space if there was none and updating the slots afterwards. */
            if ( (v_pos = get_entry(long_from_ints(xi,xj),table)) == -1) {
               v_pos = value;
               value += 9; /* we do 9 because we are adding a 
                              full 3x3 matrix                */
               add_entry(long_from_ints(xi,xj),v_pos,table);

               values_allocate_block(values,v_pos,h_pos_j,h_pos_i);
            }
            values_fill_block(values,v_pos,H_xx);
         }
         else {

            /* H_jj = (J_j)' i J_j */
            m33_A_transpose_mult_B_symmetric(temp,(*J_j),(*i));
            sm33_A_mult_B(H_jj,temp,(*J_j));

            /* b_j = (J_j)' i e_ij */
            pv_A_transpose_mult_b(b_j,temp,(*e_ij));

            /* H_ji = (J_j)' i J_i */
            m33_A_mult_B(H_xx,temp,(*J_i));

            /* H_ii = (J_i)' i J_i */
            m33_A_transpose_mult_B_symmetric(temp,(*J_i),(*i));
            sm33_A_mult_B(H_ii,temp,(*J_i));

            /* b_i = (J_i)' i e_ij */
            pv_A_transpose_mult_b(b_i,temp,(*e_ij));


            /* this block will be slightly different from the
             * previous branch, as the positions of
             * h_pos_i,h_pos_j will be flipped. We do not flip
             * xi,xj because we want the position in values of
             * competing operations to point to the same locations. */
            if ( (v_pos = get_entry(long_from_ints(xi,xj),table)) == -1) {
               v_pos = value;
               value += 9; 

               add_entry(long_from_ints(xi,xj),v_pos,table);

               values_allocate_block(values,v_pos,h_pos_i,h_pos_j);
            }
            values_fill_block(values,v_pos,H_xx);

         }

         /* insert error at position i in the vector b 
          * note that each error contains three values,
          * x, y and 0. This is the reason for the
          * multiplication by three. Additionally, 
          * because in our final system we are solving
          * for -b, we will subtract these values instead
          * of adding them                               */
         b[h_pos_i  ] -= b_i.x;
         b[h_pos_i+1] -= b_i.y;
         b[h_pos_i+2] -= b_i.t;

         /* insert error at position j in the vector b */
         b[h_pos_j  ] -= b_j.x;
         b[h_pos_j+1] -= b_j.y;
         b[h_pos_j+2] -= b_j.t;

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
            values_allocate_diagonal(values,v_pos,h_pos_i);
         }
         /* now we can insert the values. We will go by 
          * rows, so H_xx[1][2] = H_xx[4]*/
         values_fill_diagonal(values,v_pos,H_ii);

         /* repeat our actions for the matrix H_jj */
         if ( (v_pos = get_entry(long_from_ints(xj,xj),table)) == -1) {
            v_pos = value;
            value += 6; 
            add_entry(long_from_ints(xj,xj),v_pos,table);

            values_allocate_diagonal(values,v_pos,h_pos_j);
         }
         values_fill_diagonal(values,v_pos,H_jj);

      break;
      /* Our xi is a fixed node. There is no need to calculate it's
       * jacobian or the shared jacobians, as the entire row and column
       * of xi is being suppressed. Do add placeholders though.        */
      case 1:

            /* H_jj = (J_j)' i J_j */
            m33_A_transpose_mult_B_symmetric(temp,(*J_j),(*i));
            sm33_A_mult_B(H_ii,temp,(*J_j));

            /* b_i = (J_i)' i e_ij */
            pv_A_transpose_mult_b(b_j,temp,(*e_ij));

            /* insert error at position i in the vector b */
            b[h_pos_j  ] -= b_j.x;
            b[h_pos_j+1] -= b_j.y;
            b[h_pos_j+2] -= b_j.t;

            /* insert H_ii into values array */
            if ( (v_pos = get_entry(long_from_ints(xj,xj),table)) == -1) {
               v_pos = value;
               value += 6; 

               add_entry(long_from_ints(xj,xj),v_pos,table);

               values_allocate_diagonal(values,v_pos,h_pos_j);
            }
            values_fill_diagonal(values,v_pos,H_jj);

            /* insert placeholder for H_ii */
            if ( (v_pos = get_entry(long_from_ints(xi,xi),table)) == -1) {
               v_pos = value;
               value += 3; 
               add_entry(long_from_ints(xi,xi),v_pos,table);

               values_placehold_diagonal(values,v_pos,h_pos_i);
            }

      break;
      /* Our xj is a fixed node. There is no need to calculate it's
       * jacobian or the shared jacobians, as the entire row and column
       * of xj is being suppressed. Do add placeholders though.        */
      case 2:

            /* H_ii = (J_i)' i J_i */
            m33_A_transpose_mult_B_symmetric(temp,(*J_i),(*i));
            sm33_A_mult_B(H_ii,temp,(*J_i));

            /* b_i = (J_i)' i e_ij */
            pv_A_transpose_mult_b(b_i,temp,(*e_ij));

            /* insert error at position i in the vector b */
            b[h_pos_i  ] -= b_i.x;
            b[h_pos_i+1] -= b_i.y;
            b[h_pos_i+2] -= b_i.t;

            /* insert H_ii into values array */
            if ( (v_pos = get_entry(long_from_ints(xi,xi),table)) == -1) {
               v_pos = value;
               value += 6; 

               add_entry(long_from_ints(xi,xi),v_pos,table);

               values_allocate_diagonal(values,v_pos,h_pos_i);
            }
            values_fill_diagonal(values,v_pos,H_ii);

            /* insert placeholder for H_jj */
            if ( (v_pos = get_entry(long_from_ints(xj,xj),table)) == -1) {
               v_pos = value;
               value += 3; 
               add_entry(long_from_ints(xj,xj),v_pos,table);

               values_placehold_diagonal(values,v_pos,h_pos_j);
            }

      break;
      /* both of the nodes are fixed. Add placeholders and exit */
      case 3:

            /* insert placeholder for H_ii */
            if ( (v_pos = get_entry(long_from_ints(xi,xi),table)) == -1) {
               v_pos = value;
               value += 3; 
               add_entry(long_from_ints(xi,xi),v_pos,table);

               values_placehold_diagonal(values,v_pos,h_pos_i);
            }

            /* insert placeholder for H_jj */
            if ( (v_pos = get_entry(long_from_ints(xj,xj),table)) == -1) {
               v_pos = value;
               value += 3; 
               add_entry(long_from_ints(xj,xj),v_pos,table);

               values_placehold_diagonal(values,v_pos,h_pos_j);
            }

      break;
   }

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

int fix_node(pgn_handle node, pg * graph) {
   /* basically just a wrapper over a set function */
   return insert_set(node,graph->fixed_nodes);
}

void unfix_node(pgn_handle node, pg * graph) {
   /* basically just a wrapper over a set function */
   remove_set(node,graph->fixed_nodes);
}

pg * construct_pose_graph() {
   /* allocate our memory. We will start with 256
    * node and edge slots allocated              */
   pg * ret = malloc(sizeof(pg));
   ret->node = malloc(sizeof(pgn)*256);
   ret->edge = malloc(sizeof(pge)*256);
   ret->fixed_nodes = construct_set();

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
   destruct_set(graph->fixed_nodes);
   free(graph->node);
   free(graph->edge);
   free(graph);
}

int optimize(pg * graph, gsod * settings) {

   /* default number of cycles and limit
    * of granularity allowed            */
   static int default_step_limit = 0;
   static float default_cutoff = 1e-4;

   /* total number of steps allowed */
   int step_limit = default_step_limit;
   /* point at which we stop optimization */
   float cutoff = default_cutoff;

   /* the final value we return */
   int ret_val = 0;

   /* used in timing various operations */
   clock_t t_limit = CLK_MAX, t_construct = 0, t_solve = 0, t_total = clock();

   /* if we have a settings structure supplied, fill
    * out our values in preparation                 */
   if (settings != NULL) {
      /* use maximum clock value if t_limit is zero */
      t_limit = (settings->t_limit ? settings->t_limit : CLK_MAX);

      step_limit = settings->step_limit;

      /* use default cutoff if provided cutoff is zero */
      cutoff = (settings->cutoff ? settings->cutoff : default_cutoff);

      settings->steps = 0;
   }

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

   /* perform optimization loop */
   while (1) {

      if (settings != NULL) {
         settings->steps += 1;
         t_construct = clock();
      }

      bzero(b,sizeof(float)* 3 *graph->node_count);
      bzero(values,sizeof(mt)* 21 *graph->edge_count);
      clear_table(table);
      total_change = 0.0;

      int value = 0, xi, xj, b_pos;
      sm33 H_ii, H_jj;
      m33 J_i, J_j, H_ji;
      pv e_ij, b_i, b_j;
      for (int i = 0; i < graph->edge_count; ++i) {

         /* perform computation to get our jacobian blocks H_ii, H_ji, H_jj
          * and our error values b_i and b_j. These are the parameters to the
          * equations we will use the choloskey decomposition to solve       */
         compute_jacobians_and_error(&graph->edge[i],graph,&J_i,&J_j,&e_ij);

         /* compute the sections of our information matrix H that we can derive from this
          * edge, as well as the cooresponding error values of b. Add these to the
          * appropriate data structures, and return the head of the stack pointer in
          * the values array. This method takes so many arguments because this is a
          * bottleneck and one of the ways to reduce it is with some clever optimizations
          * which can only be realized by packing a lot into this one function.           */
         value = compute_constraints_and_fill_information_matrix(graph->edge[i].xi,graph->edge[i].xj,
               &J_i,&J_j,&e_ij,&graph->edge[i].information,b,values,table,graph->fixed_nodes,value);

      }

      /* fix the first entry to make our system well defined. This
       * is done by adding the identity matrix to the first diagonal */
      int fixed_val = get_entry(0,table);
      values[fixed_val  ].val += 1;
      values[fixed_val+2].val += 1;
      values[fixed_val+5].val += 1;

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

      if (settings != NULL) {
         settings->t_construct += clock() - t_construct;
         t_solve = clock();
      }

      /* perform the factorization and solve the matrix */
      perform_numerical_factorization(H,L);
      solve_system(L,b);

      if (settings != NULL)
         settings->t_solve += clock() - t_solve;

      /* collect total change */
      for (int i = 0; i < graph->node_count*3; ++i)
         total_change += fabs(b[i]);

      /* check if our decomposition has diverged. This is done to avoid letting everything
       * go to nan. Unfortunatly we cannot tell if this will happen before it does, so the
       * user may be left with poses that are partially optimized and make no sense. There
       * is not much that can be done about that, the hope is that they will go fix their
       * front end if they see this message enough. Regardless, the guarentee this section
       * provides is that a call to optimize will never result in the poses becoming undefined */
      if (isnan(total_change) || isinf(total_change) || total_change != total_change) {
         printf("[GRAPH SLAM BACKEND] optimization failure, pose graph is likely not positive semi-definite\n");
         printf("[GRAPH SLAM BACKEND] exiting prematurely, pose graph values may be trash\n");
         if (settings != NULL) {
            settings->t_total = clock() - t_total;
            settings->end_state = OES_failure;
         }
         ret_val = -1;
         goto cleanup_exit;
      }

      /* update our poses given the solution to the 
       * problem we constructed                    */
      for (int i = 0; i < graph->node_count; ++i) {
         graph->node[i].pos.x += b[i*3];
         graph->node[i].pos.y += b[(i*3)+1];
         graph->node[i].pos.t += b[(i*3)+2];
         normalize_angle(graph->node[i].pos.t);
      }

      /* check our exit conditions */
      if (clock() - t_total > t_limit) {
         if (settings != NULL) {
            settings->t_total = clock() - t_total;
            settings->end_state = OES_time_limit;
         }
         ret_val = -1;
         goto cleanup_exit;
      }

      /* if step limit was zero this will go forever */
      if (--step_limit == 0) {
         if (settings != NULL) {
            settings->t_total = clock() - t_total;
            settings->end_state = OES_step_limit;
         }
         ret_val = -1;
         goto cleanup_exit;
      }

      /* if the movement of our optimizer is 
       * below this value we assume it is done
       * and exit. Note that this may happen 
       * before the optimization converges, but
       * it will save us a few steps.          */
      if (total_change < cutoff)
         break;
   }

   /* set return values for the settings structure */
   if (settings != NULL) {
      settings->t_total = clock() - t_total;
      settings->end_state = OES_success;
   }

cleanup_exit:

   /* cleanup all the memory we allocated */
   destruct_table(table);
   free(b);
   free(values);
   free_packed_column_matrix(H);
   free_packed_column_matrix(L);

   return ret_val;

}

void print_graph_slam_optimization_data(gsod * data) {
   printf("-------------------------------------------------------------------------------\n");
   switch (data->end_state) {
      case OES_success:
         printf("result: SUCCESS\n");
         break;
      case OES_step_limit:
         printf("result: STEP LIMIT REACHED [%d]\n",
                data->step_limit);
         break;
      case OES_time_limit:
         printf("result: TIME LIMIT REACHED [%fs]\n",
                (double)data->t_limit / CLOCKS_PER_SEC);
      case OES_failure:
         printf("result: FAILURE\n");
   }
   printf("total steps: %d\n",data->steps);
   printf("total time:  %f\n",(double)data->t_total / CLOCKS_PER_SEC);
   printf("detailed time report:\n");
   printf("\t[ %14s : %20f (%%%5.2f) ]\n","construction",
         (double)data->t_construct / CLOCKS_PER_SEC, (double)data->t_construct / data->t_total);
   printf("\t[ %14s : %20f (%%%5.2f) ]\n","solving",
         (double)data->t_solve / CLOCKS_PER_SEC, (double)data->t_solve / data->t_total);
   printf("\t[ %14s : %20f (%%%5.2f) ]\n","misc",
         (double)(data->t_total - data->t_solve - data->t_construct) / CLOCKS_PER_SEC,
         (double)(data->t_total - data->t_solve - data->t_construct) / data->t_total);
   printf("-------------------------------------------------------------------------------\n");
}
