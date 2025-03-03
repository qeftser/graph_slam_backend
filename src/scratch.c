
/*
 * live testing of various
 * computations
 */

#include "graph_slam.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

int main(void) {

   /*
   float test_matrix[9] = { 4, 12, -16, 12, 37, -43, -16, -43, 98 };

   float * decomposition = basic_cholesky_factorization(test_matrix,3);

   printf("A:\n");
   print_full_matrix(test_matrix,3);
   printf("L:\n");
   print_full_matrix(decomposition,3);

   free(decomposition);
   */

   /*
   for (int i = 0; i < p->n; ++i) {
      while (row_traversal(r) != -1)
         printf("%d %d\n",r->i,r->j);
   }
   */

   /*
   int * parent = NULL;
   int nnz = construct_elimination_tree(p,&parent);
   printf("parent tree of L (%d nnz):\n",nnz);
   for (int i = 0; i < p->n; ++i)
      printf("%d ",parent[i]);
   printf("\n");

   */

   /*
   mt pcmat_test[13] = {{0,0,4.6},{2,0,1.3},{5,0,2.5},{1,1,6.4},{2,1,1.7},{5,1,3.9},{2,2,7.3},
                        {3,2,2.1},{5,2,3.1},{3,3,6.9},{4,3,2.8},{4,4,4.7},{5,5,9.9}          };
   pcmat * A = construct_packed_column_matrix(pcmat_test,6,13);

//   print_packed_column_matrix(A);

   pcmat * L = compute_symbolic_factorization(A);

   float b[6] = { 1, 1, 1, 1, 1, 1};

//   print_packed_column_matrix(L);

   perform_numerical_factorization(A,L);
   solve_system(L,b);
   printf("[%f %f %f %f %f %f ]\n",b[0],b[1],b[2],b[3],b[4],b[5]);

//   print_packed_column_matrix(L);

   free_packed_column_matrix(A);
   free_packed_column_matrix(L);
   */

   /*
   float test_matrix[9] = { 4, 12, -16, 12, 37, -43, -16, -43, 98 };
   mt test_triples[6] = {{0,0,4},{1,0,12},{1,1,37},{2,0,-16},{2,1,-43},{2,2,98}};

   A = construct_packed_column_matrix(test_triples,3,6);
   print_packed_column_matrix(A);

   L = compute_symbolic_factorization(A);

   perform_numerical_factorization(A,L);
   print_packed_column_matrix(L);

   float l[3] = { 1, 1, 1 };
   solve_system(L,l);
   printf("[ %f %f %f ]\n",l[0],l[1],l[2]);
   */

   /*
   printf("init:\n");
   pg * graph = construct_pose_graph();
   pv pos1 = { 1, 1, 0 };
   pv pos2 = { 2, 1, M_PI/2 };
   pv pos3 = { 2, 2, M_PI };
   pv pos4 = { 1, 2, -M_PI/2 };
   pv obs12 = a_from_b(pos2,pos1);
   pv obs21 = a_from_b(pos1,pos2);
   printf("observations:\n");
   printf("[%f %f %f]\n",obs12.x,obs12.y,obs12.t);
   pv obs23 = a_from_b(pos3,pos2);
   printf("[%f %f %f]\n",obs23.x,obs23.y,obs23.t);
   pv obs34 = a_from_b(pos4,pos3);
   printf("[%f %f %f]\n",obs34.x,obs34.y,obs34.t);
   pv obs41 = a_from_b(pos1,pos4);
   printf("[%f %f %f]\n",obs41.x,obs41.y,obs41.t);

   sm33 information = { 2, 1, 2, 1, 1, 2 };

   pv obs00 = { 1, 2, 3 };
   pv obs11 = { 4, 5, 6 };

   pos1.x += 10;

   printf("start:\n");
   printf("[%f %f %f]\n",pos1.x,pos1.y,pos1.t);
   printf("[%f %f %f]\n",pos2.x,pos2.y,pos2.t);
   printf("[%f %f %f]\n",pos3.x,pos3.y,pos3.t);
   printf("[%f %f %f]\n",pos4.x,pos4.y,pos4.t);
   printf("\n");

   add_node(pos1,graph);
   add_node(pos2,graph);
   add_node(pos3,graph);
   add_node(pos4,graph);

//   add_edge(0,0,&obs00,&information,graph);
//   add_edge(1,0,&obs11,&information,graph);
//   add_edge(0,1,&obs12,&information,graph);
   add_edge(0,1,&obs12,&information,graph);
   add_edge(1,2,&obs23,&information,graph);
   add_edge(2,3,&obs34,&information,graph);
   add_edge(3,0,&obs41,&information,graph);

   gsod settings;
   bzero(&settings,sizeof(gsod));

   optimize(graph,&settings);

   pos1 = graph->node[0].pos;
   pos2 = graph->node[1].pos;
   pos3 = graph->node[2].pos;
   pos4 = graph->node[3].pos;

   printf("end:\n");
   printf("[%f %f %f]\n",pos1.x,pos1.y,pos1.t);
   printf("[%f %f %f]\n",pos2.x,pos2.y,pos2.t);
   printf("[%f %f %f]\n",pos3.x,pos3.y,pos3.t);
   printf("[%f %f %f]\n",pos4.x,pos4.y,pos4.t);

   printf("\n");

   print_graph_slam_optimization_data(&settings);
   */

   pg * graph = construct_pose_graph();

// nodes: 
 pv pv0 = { 7.614383, -0.000948, 0.006676 };
 add_node(pv0,graph); 
 pv pv1 = { 7.945147, 0.000716, 0.004702 };
 add_node(pv1,graph); 
 pv pv2 = { 8.252358, 0.003277, 0.007192 };
 add_node(pv2,graph); 
 pv pv3 = { 8.558512, 0.006609, 0.010462 };
 add_node(pv3,graph); 
 pv pv4 = { 8.865058, 0.094600, 0.010454 };
 add_node(pv4,graph); 
// edges: 
// 1 -> 0	// xi: (7.945147 0.000716), 0.004702'	xj: (7.614383 -0.000948), 0.006676'
 pv   obs0 =  { -0.330769, 0.000110, 0.000585 };
 sm33 information0 = { 1.000000e+05,        
                        -0.000000e+00, 1.000000e+05,     
                        0.000000e+00, -0.000000e+00, 1.000000e+03 };
 add_edge(1,0,&obs0,&information0,graph);
// 2 -> 1	// xi: (8.252358 0.003277), 0.007192'	xj: (7.945147 0.000716), 0.004702'
 pv   obs1 =  { -0.307220, 0.000130, 0.000288 };
 sm33 information1 = { 1.000000e+05,        
                        -0.000000e+00, 1.000000e+05,     
                        0.000000e+00, -0.000000e+00, 1.000000e+03 };
 add_edge(2,1,&obs1,&information1,graph);
// 2 -> 0	// xi: (8.252358 0.003277), 0.007192'	xj: (7.614383 -0.000948), 0.006676'
 pv   obs2 =  { -0.637989, 0.000145, 0.000873 };
 sm33 information2 = { 1.000000e+05,        
                        -0.000000e+00, 1.000000e+05,     
                        0.000000e+00, -0.000000e+00, 1.000000e+03 };
 add_edge(2,0,&obs2,&information2,graph);
// 3 -> 2	// xi: (8.558512 0.006609), 0.010462'	xj: (8.252358 0.003277), 0.007192'
 pv   obs3 =  { -0.306170, 0.000134, 0.000016 };
 sm33 information3 = { 1.000000e+05,        
                        -0.000000e+00, 1.000000e+05,     
                        0.000000e+00, -0.000000e+00, 1.000000e+03 };
 add_edge(3,2,&obs3,&information3,graph);
// 3 -> 1	// xi: (8.558512 0.006609), 0.010462'	xj: (7.945147 0.000716), 0.004702'
 pv   obs4 =  { -0.613391, 0.000260, -0.017150 };
 sm33 information4 = { 1.000000e+05,        
                        -0.000000e+00, 1.000000e+05,     
                        0.000000e+00, -0.000000e+00, 1.000000e+03 };
 add_edge(3,1,&obs4,&information4,graph);
// 4 -> 3	// xi: (8.865058 0.094600), 0.010454'	xj: (8.558512 0.006609), 0.010462'
 pv   obs5 =  { -0.306975, 0.000178, 0.000008 };
 sm33 information5 = { 1.000000e+05,        
                        -0.000000e+00, 1.000000e+05,     
                        0.000000e+00, -0.000000e+00, 1.000000e+03 };
 add_edge(4,3,&obs5,&information5,graph);
// 4 -> 2	// xi: (8.865058 0.094600), 0.010454'	xj: (8.252358 0.003277), 0.007192'
 pv   obs6 =  { -0.613649, -0.088932, -0.027356 };
 sm33 information6 = { 1.000000e+05,        
                        -0.000000e+00, 1.000000e+05,     
                        0.000000e+00, -0.000000e+00, 1.000000e+03 };
 add_edge(4,2,&obs6,&information6,graph);


   gsod settings;
   bzero(&settings,sizeof(gsod));
   settings.step_limit = 5;
   int ret = optimize(graph,&settings);


   free_pose_graph(graph);

   return 0;

}
