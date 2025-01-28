
#ifndef GRAPH_SLAM

#define GRAPH_SLAM
#include "default.h"
#include "homogeneous.h"
#include "solver.h"
#include "hash.h"

/* construct a long by oring an int and a 
 * bitshifted one                         */
#define long_from_ints(i1,i2) \
   ((long)(((long)i1)<<32)|((long)i2))

/* normalize the given value as an angle
 * in radians                           */
#define normalize_angle(a) \
   ((a) = (atan2(sin((a)),cos((a)))))

/* a pointer to a node in the pose graph. Used
 * when registering new edges, so the front end
 * will want to keep track of these for it's 
 * operations                                  */
typedef unsigned int pose_graph_node_handle;
typedef pose_graph_node_handle pgn_handle;

/* a node on the pose graph. This is different
 * from the position vector to make it eaiser
 * to understand what is going on              */
typedef struct pose_graph_node {
   pv pos;
} pgn;

/* representation of a symmetric 3x3 matrix.
 * Note that there are two ways of accessing
 * the same values, one as individual elements,
 * and the other as a vector                   */
typedef struct symmetric_3x3_matrix {
   union {
      struct {
         float a00,
               a10, a11,
               a20, a21, a22;
      };
      float v[6];
   };
} sm33;

/* 3x3 matrix. Very nice. Stored
 * in row major order           */
typedef float m33[9];

/* helper method. Prints all elements
 * of a symmetric 3x3 matrix, in a 
 * decent format that is easy to read */
void print_m33(m33 * m);

/* An edge in our pose graph. Specifies
 * a constraint to use when optimizing our 
 * movement. Holds references to node in the
 * graph, as well as the actual observation
 * and it's associated information matrix. 
 * Each edge will be used to optimize the 
 * nodes. Note that the information matrix must
 * be positive semi-definite or else it will
 * cause the optimization of the pose graph
 * to diverge                                  */
typedef struct pose_graph_edge {
   pgn_handle xi;
   pgn_handle xj;
   pv observation;
   sm33 information;
} pge;

/* Data structure for holding the pose graph.
 * node and edge are simply stacks that grow
 * dynamically. These stacks are then fed into
 * the optimization function when it is called,
 * resulting in an improved set of values for
 * the nodes. It is expected that the front end
 * will access the node array pretty regularly
 * to check the current position of various nodes */
typedef struct pose_graph {
   int node_count, node_ceiling,
       edge_count, edge_ceiling;
   pgn * node;
   pge * edge;
} pg;

/* Intermediate method for computing the error function of
 * a pose graph constraint (edge) and it's resulting jacobian for
 * xi and xj.                                                     */
void compute_jacobians_and_error(pge * edge, pg * graph, m33 * J_i, m33 * J_j, pv * e_ij);
/* Intermediate method for computing sections of our information matrix and our error
 * constraint vector b given the error and jacobians we previously computed          */
void compute_constraints_and_coefficients(m33 * J_i, m33 * J_j, pv * e_ij, sm33 * information,
                                          m33 * H_ii, m33 * H_ji, m33 * H_jj, pv * b_i, pv * b_j);
/* Because we are using a sparse matrix to represent our system we need to manage a list of
 * triplets that we will shortly enter into the information matrix. This is a helper method to 
 * allow this action to not clog up our optimize method.                                       */
int insert_information_matrix_values(int xi, int xj, m33 * H_ii, m33 * H_ji, m33 * H_jj, mt * values, hash * table, int value);

/* allocate the memory and setup the
 * values for the pose graph         */
pg * construct_pose_graph();
/* deallocate all memory associated with the
 * pose graph. Note that this frees graph 
 * as well.                                 */
void free_pose_graph(pg * graph);
/* insert a node into the pose graph. This will
 * return the position in the node array that 
 * the node has been given. The user will want
 * to keep track of it in order to take advantage
 * of the optimizations performed on the node    */
pgn_handle add_node(pv pos, pg * graph);
/* insert an edge (constraint) into the pose graph. Note
 * that the information matrix must be well defined. Note 
 * also that the observation and information matricies will
 * be copied into the node, so the user can edit and resuse
 * these values without modifying the edge. This function
 * returns a handle to the observation. I am not honestly
 * sure when it would be used though. It is there if you want
 * it :)                                                      */
int        add_edge(pgn_handle xi, pgn_handle xj, pv * observation, sm33 * information, pg * graph);

/* This is the slam part of graph slam. Perform
 * optimization on the given pose graph, modifying
 * the nodes in the graph. This operation will not
 * fail if the graph is well defined, but will fail
 * if it is not. On failure, the value -1 is returned
 * ( 0 is returned on success ) and the optimization
 * is halted. The end user should check for failure, and
 * be prepared to impliment some kind of recovery method,
 * as the pose graph could have been corrupted during the
 * divergence. Note that this method may take a long time
 * to exectue if the pose graph gets large                */
int optimize(pg * graph);

#endif
