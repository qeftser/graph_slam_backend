
/* 
 * tests for the various functions
 */

#define TEST(string,check)                             \
   total++;                                            \
   if (check) {                                        \
      passing++;                                       \
      printf("%30s : \033[32mPASSED\033[0m\n",string); \
   }                                                   \
   else                                                \
      printf("%30s : \033[31mFAILED\033[0m\n",string)

#include "homogeneous.h"
#include "solver.h"
#include "hash.h"
#include "graph_slam.h"
#include "set.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>

/* helper method for computing the exact observation
 * given from seeing value a from value b. Used to
 * generate parameters for passing to our pose graph
 * during testing                                    */
pv a_from_b(pv a, pv b) {
   ht hta = as_homogeneous_transformation(a);
   ht htb = as_homogeneous_transformation(b);
   ht htb_inv = invert_homogeneous_transformation(&htb);
   ht a_f_b = merge_homogeneous_transformation(&htb_inv,&hta);
   return destruct_homogeneous_transformation(&a_f_b);
}

/* a function for checking the equality of a list of 
 * double values. Count is the number of pairs expected,
 * and epsilon is the tolarated error in the matches. 
 * After these two values, (count) pairs of doubles are 
 * expected to be passed to the function, resulting in
 * 2*(count)+2 total arguments. A 1 is returned if all
 * values are seen as equal, and a 0 is returned otherwise */
int check_equality_double(double epsilon, int count, ...) {
   va_list args;
   va_start(args,count);
   int res = 1;
   for (int i = 0; i < count; ++i) {
      double dbl_1 = va_arg(args,double);
      double dbl_2 = va_arg(args,double);
      if (fabs(dbl_1 - dbl_2) > epsilon)
         res = 0;
   }
   va_end(args);
   return res;
}

int main(void) {

   int passing = 0;
   int total   = 0;

   printf("homogeneous.h tests\n");
   {
      ht t1;
      t1.c1 = 1; t1.s1 = 2; t1.tx = 3;
      t1.s2 = 4; t1.c2 = 5; t1.ty = 6;
      t1.z1 = 7; t1.z2 = 8; t1.w  = 9;
      ht t2;
      t2.v[0] = 1; t2.v[1] = 2; t2.v[2] = 3;
      t2.v[3] = 4; t2.v[4] = 5; t2.v[5] = 6;
      t2.v[6] = 7; t2.v[7] = 8; t2.v[8] = 9;
      TEST("memory layout check",(0 == memcmp(&t1,&t2,sizeof(ht))));
   }
   {
      hv v1 = { 1, 2, 1 };
      hv v2 = construct_homogeneous_vector(1,2);
      TEST("vector construction test",(0 == memcmp(&v1,&v2,sizeof(hv))));
   }
   {
      pv pos = { 1, 2, 0.5 };
      ht t1 = as_homogeneous_transformation(pos);
      pv rpos = destruct_homogeneous_transformation(&t1);
      TEST("t2v and v2t tests",(0 == memcmp(&pos,&rpos,sizeof(pv))));
   }
   {
      ht t1 = construct_homogeneous_transformation(0.5,1,2);
      ht t2 = construct_homogeneous_transformation(0.25,3,4);
      ht t3 = invert_homogeneous_transformation(&t1);
      ht t4 = merge_homogeneous_transformation(&t1,&t2);
      ht t5 = merge_homogeneous_transformation(&t3,&t4);
      normalize_homogeneous_transformation(&t5);
      TEST("merge and inversion checks",check_equality_double(0.0001,9,t2.v[0],t5.v[0],
               t2.v[1],t5.v[1],t2.v[2],t5.v[2],t2.v[3],t5.v[3],t2.v[4],t5.v[4],t2.v[5],
               t5.v[5],t2.v[6],t5.v[6],t2.v[7],t5.v[7],t2.v[8],t5.v[8]));
   }
   {
      ht t1 = construct_homogeneous_transformation(30,1,2);
      ht t2 = invert_homogeneous_transformation(&t1);
      ht t3 = invert_homogeneous_transformation(&t2);
      TEST("double inversion",(fabs(t1.v[0] - t3.v[0]) < 0.0001) &&
                              (fabs(t1.v[1] - t3.v[1]) < 0.0001) &&
                              (fabs(t1.v[2] - t3.v[2]) < 0.0001) &&
                              (fabs(t1.v[3] - t3.v[3]) < 0.0001) &&
                              (fabs(t1.v[4] - t3.v[4]) < 0.0001) &&
                              (fabs(t1.v[5] - t3.v[5]) < 0.0001) &&
                              (fabs(t1.v[6] - t3.v[6]) < 0.0001) &&
                              (fabs(t1.v[7] - t3.v[7]) < 0.0001) &&
                              (fabs(t1.v[8] - t3.v[8]) < 0.0001));
   }

   printf("solver.h tests\n");
   {
      mt triplets[13] = {{0,0,4.6},{2,0,1.3},{5,0,2.5},{1,1,6.4},{2,1,1.7},{5,1,3.9},{2,2,7.3},
                        {3,2,2.1},{5,2,3.1},{3,3,6.9},{4,3,2.8},{4,4,4.7},{5,5,9.9}          };
      int expected_colp[7] = { 0, 3, 6, 9, 11, 12, 13 };
      int expected_rx[13] = { 0, 2, 5, 1, 2, 5, 2, 3, 5, 3, 4, 4, 5 };
      float expected_val[13] = { 4.6, 1.3, 2.5, 6.4, 1.7, 3.9, 7.3, 2.1, 3.1, 6.9, 2.8, 4.7, 9.9 };
      pcmat * A = construct_packed_column_matrix(triplets,6,13);
      TEST("pcmat construction",(A->n == 6)                                        &&
                                (A->nnz == 13)                                     &&
                                (memcmp(A->colp,expected_colp,sizeof(int)*7) == 0) &&
                                (memcmp(A->rx,expected_rx,sizeof(int)*13) == 0)    &&
                                (memcmp(A->val,expected_val,sizeof(float)*13) == 0));
      free_packed_column_matrix(A);
   }
   {
      float mat[9] = {  4,  12, -16,
                       12,  37, -43,
                      -16, -43,  98 };
      mt triplets[6] = {{0,0,4},{1,0,12},{1,1,37},{2,0,-16},{2,1,-43},{2,2,98}};

      pcmat * A = construct_packed_column_matrix(triplets,3,6);
      pcmat * L = compute_symbolic_factorization(A);
      perform_numerical_factorization(A,L);

      float b[3] = { 21, 63, 189 };
      solve_system(L,b);
      float res[3] = { mat[0] * b[0] + mat[1] * b[1] + mat[2] * b[2],
                       mat[3] * b[0] + mat[4] * b[1] + mat[5] * b[2],
                       mat[6] * b[0] + mat[7] * b[1] + mat[8] * b[2] };
      TEST("dense decomposition",(fabs(21.0  - res[0]) < 0.001) &&
                                 (fabs(63.0  - res[1]) < 0.001) &&
                                 (fabs(189.0 - res[2]) < 0.001));
      free_packed_column_matrix(A);
      free_packed_column_matrix(L);
   }
   {
      float mat[36] = { 4.6, 0.0, 1.3, 0.0, 0.0, 2.5,
                        0.0, 6.4, 1.7, 0.0, 0.0, 3.9,
                        1.3, 1.7, 7.3, 2.1, 0.0, 3.1,
                        0.0, 0.0, 2.1, 6.9, 2.8, 0.0,
                        0.0, 0.0, 0.0, 2.8, 4.7, 0.0,
                        2.5, 3.9, 3.1, 0.0, 0.0, 9.9 };
      mt triplets[13] = {{0,0,4.6},{2,0,1.3},{5,0,2.5},{1,1,6.4},{2,1,1.7},{5,1,3.9},{2,2,7.3},
                        {3,2,2.1},{5,2,3.1},{3,3,6.9},{4,3,2.8},{4,4,4.7},{5,5,9.9}          };

      pcmat * A = construct_packed_column_matrix(triplets,6,13);
      pcmat * L = compute_symbolic_factorization(A);
      perform_numerical_factorization(A,L);

      float b[6] = { 0.5, 1.5, 2.75, 3.5, 4.25, 5.5 };
      solve_system(L,b);
      float res[6] = { mat[ 0]*b[0] + mat[ 1]*b[1] + mat[ 2]*b[2] + mat[ 3]*b[3] + mat[ 4]*b[4] + mat[ 5]*b[5],
                       mat[ 6]*b[0] + mat[ 7]*b[1] + mat[ 8]*b[2] + mat[ 9]*b[3] + mat[10]*b[4] + mat[11]*b[5],
                       mat[12]*b[0] + mat[13]*b[1] + mat[14]*b[2] + mat[15]*b[3] + mat[16]*b[4] + mat[17]*b[5],
                       mat[18]*b[0] + mat[19]*b[1] + mat[20]*b[2] + mat[21]*b[3] + mat[22]*b[4] + mat[23]*b[5],
                       mat[24]*b[0] + mat[25]*b[1] + mat[26]*b[2] + mat[27]*b[3] + mat[28]*b[4] + mat[29]*b[5],
                       mat[30]*b[0] + mat[31]*b[1] + mat[32]*b[2] + mat[33]*b[3] + mat[34]*b[4] + mat[35]*b[5] };
      TEST("sparse decomposition",(fabs(0.5  - res[0]) < 0.001) &&
                                  (fabs(1.5  - res[1]) < 0.001) &&
                                  (fabs(2.75 - res[2]) < 0.001) &&
                                  (fabs(3.5  - res[3]) < 0.001) &&
                                  (fabs(4.25 - res[4]) < 0.001) &&
                                  (fabs(5.5  - res[5]) < 0.001));
      free_packed_column_matrix(A);
      free_packed_column_matrix(L);
   }
   printf("hash.h tests\n");
   {
      hash * table = construct_table(10);
      TEST("table construct/destruct",(table->table_size == 10) &&
                                      (table->entry_ceiling == 10) &&
                                      (table->entry == 0));
      destruct_table(table);
   }
   {
      hash * table = construct_table(20);
      add_entry(1,2,table);
      add_entry(3,4,table);
      add_entry(5,6,table);
      TEST("table add/get",(get_entry(1,table) == 2) &&
                           (get_entry(3,table) == 4) &&
                           (get_entry(5,table) == 6));
      destruct_table(table);
   }
   {
      hash * table = construct_table(5);
      for (int i = 0; i < 30; ++i) {
         add_entry(i,30-i,table);
      }
      TEST("grow past initial size",(get_entry(15,table) == 15) &&
                                    (get_entry(10,table) == 20) &&
                                    (get_entry(20,table) == 10));
      destruct_table(table);
   }
   {
      hash * table = construct_table(1); /* try and ensure we will traverse a chain */
      add_entry(1,2,table); 
      add_entry(3,4,table);
      add_entry(5,6,table);
      TEST("table get non-entry 1",(get_entry(7,table) == -1));
      destruct_table(table);
   }
   {
      hash * table = construct_table(1); 
      TEST("table get non-entry 1",(get_entry(5,table) == -1));
      destruct_table(table);
   }
   {
      hash * table = construct_table(10);
      add_entry(1,2,table); 
      add_entry(3,4,table);
      add_entry(5,6,table);
      clear_table(table);
      TEST("table clear",(get_entry(1,table) == -1) &&
                         (get_entry(3,table) == -1) &&
                         (get_entry(5,table) == -1));
      destruct_table(table);
   }
   printf("set.h tests\n");
   {
      set * s = construct_set();
      destruct_set(s);
      TEST("construct/destruct",1);
   }
   {
      set * s = construct_set();
      insert_set(5,s);
      TEST("insert 1",member_set(5,s) &&
                     !member_set(6,s));
      destruct_set(s);
   }
   {
      set * s = construct_set();
      insert_set(1,s);
      insert_set(5,s);
      insert_set(12,s);
      insert_set(23,s);
      insert_set(100,s);
      insert_set(84,s);
      insert_set(65,s);
      TEST("insert 2",member_set(1,s) &&
                      member_set(5,s) &&
                      member_set(12,s) &&
                      member_set(23,s) &&
                      member_set(100,s) &&
                      member_set(84,s) &&
                      member_set(65,s) &&
                     !member_set(6,s));
      destruct_set(s);
   }
   {
      set * s = construct_set();
      insert_set(50000,s);
      TEST("insert expand",member_set(50000,s) &&
                          !member_set(50001,s) &&
                          !member_set(49999,s));
      destruct_set(s);
   }
   {
      set * s = construct_set();
      insert_set(1,s);
      insert_set(5,s);
      insert_set(12,s);
      insert_set(23,s);
      insert_set(100,s);
      insert_set(84,s);
      insert_set(65,s);
      remove_set(1,s);
      remove_set(5,s);
      remove_set(12,s);
      remove_set(23,s);
      remove_set(100,s);
      remove_set(84,s);
      remove_set(65,s);
      TEST("remove",!member_set(1,s) &&
                    !member_set(5,s) &&
                    !member_set(12,s) &&
                    !member_set(23,s) &&
                    !member_set(100,s) &&
                    !member_set(84,s) &&
                    !member_set(65,s) &&
                    !member_set(6,s));
      destruct_set(s);
   }
   printf("graph_slam.h tests\n");
   {
      pg * graph = construct_pose_graph();
      free_pose_graph(graph);
      TEST("construct/destruct",1);
   }
   {
      pg * graph = construct_pose_graph();
      pv pos1 = { 1, 1, 0 };
      pv pos2 = { 2, 1, M_PI/2 };
      pv pos3 = { 2, 2, M_PI };
      pv pos4 = { 1, 2, -M_PI/2 };
      pv obs12 = a_from_b(pos2,pos1);
      pv obs23 = a_from_b(pos3,pos2);
      pv obs34 = a_from_b(pos4,pos3);
      pv obs41 = a_from_b(pos1,pos4);

      sm33 information = { 2, -1, 2, -1, -1, 2 };

      add_node(pos1,graph);
      add_node(pos2,graph);
      add_node(pos3,graph);
      add_node(pos4,graph);

      add_edge(0,1,&obs12,&information,graph);
      add_edge(1,2,&obs23,&information,graph);
      add_edge(2,3,&obs34,&information,graph);
      add_edge(3,0,&obs41,&information,graph);

      TEST("add nodes+edges",(graph->node_count == 4) &&
                             (graph->edge_count == 4) &&
                             (memcmp(&graph->node[0].pos,&pos1,sizeof(pv)) == 0) &&
                             (memcmp(&graph->node[1].pos,&pos2,sizeof(pv)) == 0) &&
                             (memcmp(&graph->node[2].pos,&pos3,sizeof(pv)) == 0) &&
                             (memcmp(&graph->node[3].pos,&pos4,sizeof(pv)) == 0) &&
                             (graph->edge[0].xi == 0) && (graph->edge[0].xj == 1) &&
                             (graph->edge[1].xi == 1) && (graph->edge[1].xj == 2) &&
                             (graph->edge[2].xi == 2) && (graph->edge[2].xj == 3) &&
                             (graph->edge[3].xi == 3) && (graph->edge[3].xj == 0) &&
                             (memcmp(&graph->edge[0].observation,&obs12,sizeof(pv)) == 0) &&
                             (memcmp(&graph->edge[1].observation,&obs23,sizeof(pv)) == 0) &&
                             (memcmp(&graph->edge[2].observation,&obs34,sizeof(pv)) == 0) &&
                             (memcmp(&graph->edge[3].observation,&obs41,sizeof(pv)) == 0) &&
                             (memcmp(&graph->edge[0].information,&information,sizeof(sm33)) == 0) &&
                             (memcmp(&graph->edge[1].information,&information,sizeof(sm33)) == 0) &&
                             (memcmp(&graph->edge[2].information,&information,sizeof(sm33)) == 0) &&
                             (memcmp(&graph->edge[3].information,&information,sizeof(sm33)) == 0));
      free_pose_graph(graph);
   }
   {
      pg * graph = construct_pose_graph();
      pv pos1 = { 1, 1, 0 };
      pv pos2 = { 2, 1, M_PI/2 };
      pv obs12 = a_from_b(pos2,pos1);

      sm33 information = { 2, -1, 2, -1, -1, 2 };

      add_node(pos1,graph);
      add_node(pos2,graph);
      add_edge(0,1,&obs12,&information,graph);

      m33 J_i, J_j;
      pv e_ij;
      compute_jacobians_and_error(&graph->edge[0],graph,&J_i,&J_j,&e_ij);

      TEST("ensure jacobian+error",(fabs(e_ij.x) < 0.00001) && (fabs(e_ij.y) < 0.00001) &&
                                   (fabs(e_ij.t) < 0.00001) && (fabs(J_i[0]) < 0.00001) &&
                                   (fabs(J_i[2]) < 0.00001) && (fabs(J_i[4]) < 0.00001) &&
                                   (fabs(J_i[6]) < 0.00001) && (fabs(J_i[7]) < 0.00001) &&
                                   (fabs(J_j[0]) < 0.00001) && (fabs(J_j[2]) < 0.00001) &&
                                   (fabs(J_j[4]) < 0.00001) && (fabs(J_j[5]) < 0.00001) &&
                                   (fabs(J_j[6]) < 0.00001) && (fabs(J_j[7]) < 0.00001) &&
                                   (fabs(-1 - J_i[1]) < 0.00001) && (fabs(1 - J_i[3]) < 0.00001) &&
                                   (fabs(-1 - J_i[5]) < 0.00001) && (fabs(-1 - J_i[8]) < 0.00001) &&
                                   (fabs(1 - J_j[1]) < 0.00001) && (fabs(-1 - J_j[3]) < 0.00001) &&
                                   (fabs(1 - J_j[8]) < 0.00001));
      free_pose_graph(graph);
   }
   {
      pg * graph = construct_pose_graph();
      pv pos1 = { 1, 1, 0 };
      pv pos2 = { 2, 1, M_PI/2 };
      pv obs12 = a_from_b(pos2,pos1);

      sm33 information = { 2, -1, 2, -1, -1, 2 };

      add_node(pos1,graph);
      add_node(pos2,graph);
      add_edge(0,1,&obs12,&information,graph);

      int b_pos, v_pos, xi, xj;
      m33 J_i, J_j, H_ii, H_ji, H_jj;
      pv e_ij, b_i, b_j;
      compute_jacobians_and_error(&graph->edge[0],graph,&J_i,&J_j,&e_ij);
      compute_constraints_and_coefficients(&J_i,&J_j,&e_ij,&graph->edge[0].information,
                                           &H_ii,&H_ji,&H_jj,&b_i,&b_j);
      TEST("ensure H and b values",(fabs(H_ii[0] - 2) < 0.00001) && (fabs(H_ii[1] - 1) < 0.00001) &&
                                   (fabs(H_ii[2] + 1) < 0.00001) && (fabs(H_ii[3] - 1) < 0.00001) &&
                                   (fabs(H_ii[4] - 2) < 0.00001) && (fabs(H_ii[5] + 2) < 0.00001) &&
                                   (fabs(H_ii[6] + 1) < 0.00001) && (fabs(H_ii[7] + 2) < 0.00001) &&
                                   (fabs(H_ii[8] - 2) < 0.00001) && (fabs(H_jj[0] - 2) < 0.00001) &&
                                   (fabs(H_jj[1] - 1) < 0.00001) && (fabs(H_jj[2] - 1) < 0.00001) &&
                                   (fabs(H_jj[3] - 1) < 0.00001) && (fabs(H_jj[4] - 2) < 0.00001) &&
                                   (fabs(H_jj[5] + 1) < 0.00001) && (fabs(H_jj[6] - 1) < 0.00001) &&
                                   (fabs(H_jj[7] + 1) < 0.00001) && (fabs(H_jj[8] - 2) < 0.00001) &&
                                   (fabs(H_ji[0] + 2) < 0.00001) && (fabs(H_ji[1] + 1) < 0.00001) &&
                                   (fabs(H_ji[2] - 1) < 0.00001) && (fabs(H_ji[3] + 1) < 0.00001) &&
                                   (fabs(H_ji[4] + 2) < 0.00001) && (fabs(H_ji[5] - 2) < 0.00001) &&
                                   (fabs(H_ji[6] + 1) < 0.00001) && (fabs(H_ji[7] - 1) < 0.00001) &&
                                   (fabs(H_ji[8] + 1) < 0.00001) && (fabs(b_i.x) < 0.00001) && 
                                   (fabs(b_i.y) < 0.00001) && (fabs(b_i.t) < 0.00001) &&
                                   (fabs(b_j.x) < 0.00001) && (fabs(b_j.y) < 0.00001) && (fabs(b_j.t) < 0.00001));
      free_pose_graph(graph);
   }
   {
      pg * graph = construct_pose_graph();
      pv pos1 = { 1, 1, 0 };
      pv pos2 = { 2, 1, M_PI/2 };
      pv obs12 = a_from_b(pos2,pos1);

      sm33 information = { 2, -1, 2, -1, -1, 2 };

      add_node(pos1,graph);
      add_node(pos2,graph);
      add_edge(0,1,&obs12,&information,graph);

      int value = 0, b_pos, v_pos;
      m33 J_i, J_j, H_ii, H_ji, H_jj;
      pv e_ij, b_i, b_j;
      compute_jacobians_and_error(&graph->edge[0],graph,&J_i,&J_j,&e_ij);
      compute_constraints_and_coefficients(&J_i,&J_j,&e_ij,&graph->edge[0].information,
                                           &H_ii,&H_ji,&H_jj,&b_i,&b_j);

      hash  * table = construct_table(graph->edge_count);
      mt * values = malloc(sizeof(mt)* 21 *graph->edge_count);
      bzero(values,sizeof(mt)* 21 *graph->edge_count);
      clear_table(table);

      value = insert_information_matrix_values(0,1,&H_ii,&H_ji,&H_jj,values,table,value);
      pcmat * H = construct_packed_column_matrix(values,graph->node_count*3,value);
      TEST("ensure correct H structure",(value == 21) &&
                                        (fabs(H->val[0] - 2) < 0.00001) && (fabs(H->val[1] - 1) < 0.00001) &&
                                        (fabs(H->val[2] + 1) < 0.00001) && (fabs(H->val[3] + 2) < 0.00001) &&
                                        (fabs(H->val[4] + 1) < 0.00001) && (fabs(H->val[5] + 1) < 0.00001) &&
                                        (fabs(H->val[6] - 2) < 0.00001) && (fabs(H->val[7] + 2) < 0.00001) &&
                                        (fabs(H->val[8] + 1) < 0.00001) && (fabs(H->val[9] + 2) < 0.00001) &&
                                        (fabs(H->val[10] - 1) < 0.00001) && (fabs(H->val[11] - 2) < 0.00001) &&
                                        (fabs(H->val[12] - 1) < 0.00001) && (fabs(H->val[13] - 2) < 0.00001) &&
                                        (fabs(H->val[14] + 1) < 0.00001) && (fabs(H->val[15] - 2) < 0.00001) &&
                                        (fabs(H->val[16] - 1) < 0.00001) && (fabs(H->val[17] - 1) < 0.00001) &&
                                        (fabs(H->val[18] - 2) < 0.00001) && (fabs(H->val[19] + 1) < 0.00001) &&
                                        (fabs(H->val[20] - 2) < 0.00001));
      destruct_table(table);
      free_packed_column_matrix(H);
      free(values);
      free_pose_graph(graph);
   }
   {
      pg * graph = construct_pose_graph();
      pv poses[4] = {{1, 1, 0},{ 2, 1, M_PI/2},{2, 2, M_PI},{1, 2, -M_PI/2}};
      pv obs[4] = {a_from_b(poses[1],poses[0]),a_from_b(poses[2],poses[1]),
                   a_from_b(poses[3],poses[2]),a_from_b(poses[0],poses[3])};
      sm33 information = { 2, 1, 2, 1, 1, 2 };

      poses[3].x += 10;
      
      for (int i = 0; i < 4; ++i) {
         add_node(poses[i],graph);
         add_edge(i,(i+1)%4,&obs[i],&information,graph);
      }

      optimize(graph,NULL);

      for(int i = 0; i < 4; ++i)
         poses[i] = graph->node[i].pos;

      TEST("optimize test 1",check_equality_double(0.0001,12,poses[0].x,1.0,poses[0].y,1.0,poses[0].t,0.0,
                                                             poses[1].x,2.0,poses[1].y,1.0,poses[1].t,M_PI/2,
                                                             poses[2].x,2.0,poses[2].y,2.0,poses[2].t,M_PI,
                                                             poses[3].x,1.0,poses[3].y,2.0,poses[3].t,-M_PI/2));
      free_pose_graph(graph);
   }
   {
      pg * graph = construct_pose_graph();
      pv poses[4] = {{1, 1, 0},{ 2, 1, M_PI/2},{2, 2, M_PI},{1, 2, -M_PI/2}};
      pv obs[4] = {a_from_b(poses[1],poses[0]),a_from_b(poses[2],poses[1]),
                   a_from_b(poses[3],poses[2]),a_from_b(poses[0],poses[3])};
      sm33 information = { 2, 1, 2, 1, 1, 2 };
      
      poses[0].x += 10;

      for (int i = 0; i < 4; ++i) {
         add_node(poses[i],graph);
         add_edge(i,(i+1)%4,&obs[i],&information,graph);
      }

      optimize(graph,NULL);

      for(int i = 0; i < 4; ++i)
         poses[i] = graph->node[i].pos;

      TEST("optimize test 2",check_equality_double(0.0001,12,poses[0].x,11.0,poses[0].y,1.0,poses[0].t,0.0,
                                                             poses[1].x,12.0,poses[1].y,1.0,poses[1].t,M_PI/2,
                                                             poses[2].x,12.0,poses[2].y,2.0,poses[2].t,M_PI,
                                                             poses[3].x,11.0,poses[3].y,2.0,poses[3].t,-M_PI/2));
      free_pose_graph(graph);
   }


   printf("total: %d tests, %d passed (%.2f%%)\n",total,passing,((double)passing/total)*100);

   return 0;
}
