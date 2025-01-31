
#ifndef GRAPH_SLAM

#define GRAPH_SLAM
#include "default.h"
#include "homogeneous.h"
#include "solver.h"
#include "hash.h"
#include "set.h"
#include <time.h>

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

/* multiply the m33 A transpose by the
 * sm33 symmetric matrix B and store
 * the result in res, transpose left during
 * multiplication. This is a helper method
 * for compute_constraints_and_fill_information_matrix */
#define m33_A_transpose_mult_B_symmetric(res_m33,A,B) \
{                                                     \
   res_m33[0] = A[0]*B.a00 + A[3]*B.a10 + A[6]*B.a20; \
   res_m33[1] = A[0]*B.a10 + A[3]*B.a11 + A[6]*B.a21; \
   res_m33[2] = A[0]*B.a20 + A[3]*B.a21 + A[6]*B.a22; \
   res_m33[3] = A[1]*B.a00 + A[4]*B.a10 + A[7]*B.a20; \
   res_m33[4] = A[1]*B.a10 + A[4]*B.a11 + A[7]*B.a21; \
   res_m33[5] = A[1]*B.a20 + A[4]*B.a21 + A[7]*B.a22; \
   res_m33[6] = A[2]*B.a00 + A[5]*B.a10 + A[8]*B.a20; \
   res_m33[7] = A[2]*B.a10 + A[5]*B.a11 + A[8]*B.a21; \
   res_m33[8] = A[2]*B.a20 + A[5]*B.a21 + A[8]*B.a22; \
} 

/* multiply two m33 matricies A (left side)
 * and B (right side), storing the result in
 * a sm33 symmetric matrix. This assumes that
 * the multiplication A * B will produce a
 * symmetric matrix. This is a helper method
 * for comoute_constraints_and_fill_information_matrix */
#define sm33_A_mult_B(res_sm33,A,B)                  \
{                                                    \
   res_sm33.a00 = A[0]*B[0] + A[1]*B[3] + A[2]*B[6]; \
   res_sm33.a10 = A[3]*B[0] + A[4]*B[3] + A[5]*B[6]; \
   res_sm33.a11 = A[3]*B[1] + A[4]*B[4] + A[5]*B[7]; \
   res_sm33.a20 = A[6]*B[0] + A[7]*B[3] + A[8]*B[6]; \
   res_sm33.a21 = A[6]*B[1] + A[7]*B[4] + A[8]*B[7]; \
   res_sm33.a22 = A[6]*B[2] + A[7]*B[5] + A[8]*B[8]; \
}

#define m33_A_mult_B(res_m33,A,B)                  \
{                                                  \
   res_m33[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6]; \
   res_m33[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7]; \
   res_m33[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8]; \
   res_m33[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6]; \
   res_m33[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7]; \
   res_m33[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8]; \
   res_m33[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6]; \
   res_m33[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7]; \
   res_m33[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8]; \
}

/* multiply the m33 A transpose by the
 * pv vector b, producing a new pv. This is
 * is a helper method for 
 * compute_constraints_and_fill_information_matrix */
#define pv_A_transpose_mult_b(res_pv,A,b)     \
{                                             \
   res_pv.x = A[0]*b.x + A[1]*b.y + A[2]*b.t; \
   res_pv.y = A[3]*b.x + A[4]*b.y + A[5]*b.t; \
   res_pv.t = A[6]*b.x + A[7]*b.y + A[8]*b.t; \
}

/* setup space in an array of mt matrix_triplets 
 * (values) by setting the diagonal sm33 values 
 * at v_pos to the offset of pos. This is a partial
 * matinence function used in
 * compute_constraints_and_fill_information_matrix */
#define values_allocate_diagonal(values,v_pos,pos)           \
{                                                            \
   values[v_pos  ].col = pos  ; values[v_pos  ].row = pos  ; \
   values[v_pos+1].col = pos  ; values[v_pos+1].row = pos+1; \
   values[v_pos+2].col = pos+1; values[v_pos+2].row = pos+1; \
   values[v_pos+3].col = pos  ; values[v_pos+3].row = pos+2; \
   values[v_pos+4].col = pos+1; values[v_pos+4].row = pos+2; \
   values[v_pos+5].col = pos+2; values[v_pos+5].row = pos+2; \
}

/* setup and fill space in any array of mt matrix_triplets
 * with values of 0 on pos, pos+1, and pos+2 positions on
 * the diagonal. This method is a helper used to define fixed
 * nodes in our information matrix. Because the matrix  is sparse,
 * we must have values on the diagonal for the solver to function
 * properly. This method does the work in one line in our 
 * functions, which makes them eaiser to comprehend              */
#define values_placehold_diagonal(values,v_pos,pos)   \
{                                                     \
   values[v_pos  ].col = values[v_pos  ].row = pos  ; \
   values[v_pos+1].col = values[v_pos+1].row = pos+1; \
   values[v_pos+2].col = values[v_pos+2].row = pos+2; \
   values[v_pos  ].val = 0;                           \
   values[v_pos+1].val = 0;                           \
   values[v_pos+2].val = 0;                           \
}

/* setup spance in an array of mt matrix_triplets (values)
 * by setting the row and column parameters of a 3x3 
 * block (9 elements). We will fill these values later.
 * We need two positions because we are not on the diagonal,
 * we are in an arbitrary position below it.                */
#define values_allocate_block(values,v_pos,pos_i,pos_j)              \
{                                                                    \
   values[v_pos  ].col = pos_i  ; values[v_pos  ].row = pos_j  ; \
   values[v_pos+1].col = pos_i+1; values[v_pos+1].row = pos_j  ; \
   values[v_pos+2].col = pos_i+2; values[v_pos+2].row = pos_j  ; \
   values[v_pos+3].col = pos_i  ; values[v_pos+3].row = pos_j+1; \
   values[v_pos+4].col = pos_i+1; values[v_pos+4].row = pos_j+1; \
   values[v_pos+5].col = pos_i+2; values[v_pos+5].row = pos_j+1; \
   values[v_pos+6].col = pos_i  ; values[v_pos+6].row = pos_j+2; \
   values[v_pos+7].col = pos_i+1; values[v_pos+7].row = pos_j+2; \
   values[v_pos+8].col = pos_i+2; values[v_pos+8].row = pos_j+2; \
}

/* Add the values of the symmetric matrix A_sm33 to the
 * diagonal entries at v_pos in values. This macro should
 * be called in combination with values_allocate_diagonal.
 * That macro allocates the space, and this one fills it  */
#define values_fill_diagonal(values,v_pos,A_sm33) \
{                                                 \
   values[v_pos  ].val += A_sm33.a00;             \
   values[v_pos+1].val += A_sm33.a10;             \
   values[v_pos+2].val += A_sm33.a11;             \
   values[v_pos+3].val += A_sm33.a20;             \
   values[v_pos+4].val += A_sm33.a21;             \
   values[v_pos+5].val += A_sm33.a22;             \
}

/* Add the values of the 3x3 matrix A_m33 to the block held
 * in the values array at v_pos. This macro is meant to be
 * called in combination with values_allocate_diagonal. This
 * is the one that fills val, as opposed to col and row.    */
#define values_fill_block(values,v_pos,A_m33) \
{                                             \
   values[v_pos  ].val += A_m33[0];           \
   values[v_pos+1].val += A_m33[1];           \
   values[v_pos+2].val += A_m33[2];           \
   values[v_pos+3].val += A_m33[3];           \
   values[v_pos+4].val += A_m33[4];           \
   values[v_pos+5].val += A_m33[5];           \
   values[v_pos+6].val += A_m33[6];           \
   values[v_pos+7].val += A_m33[7];           \
   values[v_pos+8].val += A_m33[8];           \
}

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
 * to check the current position of various nodes.
 * The set fixed_nodes is used to keep track of 
 * nodes that are fixed and will not be adjusted
 * in optimization.                             */
typedef struct pose_graph {
   int node_count, node_ceiling,
       edge_count, edge_ceiling;
   pgn * node;
   pge * edge;
   set * fixed_nodes;
} pg;

/* Intermediate method for computing the error function of
 * a pose graph constraint (edge) and it's resulting jacobian for
 * xi and xj.                                                     */
void compute_jacobians_and_error(pge * edge, pg * graph, m33 * J_i, m33 * J_j, pv * e_ij);

/* Intermediate method for computing sections of our information matrix and our error
 * constraint vector b given the error and jacobians we previously computed          */
/* Because we are using a sparse matrix to represent our system we need to manage a list of
 * triplets that we will shortly enter into the information matrix. This is a helper method to 
 * allow this action to not clog up our optimize method.                                       */
/* Does all the math to combine the jacobians, error function, and information matrix into a 
 * system with H as the overall information matrix and b as the constraint vector.           */
int compute_constraints_and_fill_information_matrix(int xi, int xj, m33 * J_i, m33 * J_j, pv * e_ij, sm33 * i,
                                                    float * b, mt * values, hash * table, set * fixed,
                                                    int value);
/* allocate the memory and setup the
 * values for the pose graph         */
pg * construct_pose_graph(void);
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
 * it :)                                                      
 * Note: this function returns -1 if a memory allocation
 * error occured                                             */
int        add_edge(pgn_handle xi, pgn_handle xj, pv * observation, sm33 * information, pg * graph);
/* mark a node in the pose graph as fixed. This means that 
 * the node will not be changed in the optimization process.
 * Note: this function returns -1 if a memory allocation 
 * error occured, and zero otherwise. */
int fix_node(pgn_handle node, pg * graph);
/* remove the fixed status on a node. This will allow the 
 * node to be modified in the optimize function          */
void unfix_node(pgn_handle node, pg * graph);

/* maximum value the clock_t variable will hold */
#define CLK_MAX (~(((clock_t)1)<<((((sizeof(clock_t))*8))-1)))

/* Set of values the for the end state of the optimization
 * function. succcess is set if the optimization was allowed 
 * to finish normally, failure is set on divergence, step_limit
 * is set if the step limit was hit before optimization finished,
 * and time_limit is set if the time limit was hit before
 * optimization finished                                         */
enum optimize_end_state{ OES_success,  OES_step_limit,  OES_time_limit,  OES_failure};

/* parameters and fields to pass to the 
 * optimization function to control it's
 * behavior and to get more data out of it.
 * Otherwise the solver is run with it's
 * default parameters and settings. Values 
 * labeled as input are set by the user and 
 * read by the optimize function, and values
 * labeled as output are set by the optimize 
 * function to be read by the user           */
typedef struct graph_slam_optimization_data {
   clock_t t_total;     /*OUTPUT how many clocks did the optimization run for?              */
   clock_t t_construct; /*OUTPUT how many clocks were spent on construction of the problem? */
   clock_t t_solve;     /*OUTPUT how many clock cycles were spent on solving the system?    */
   int     steps;       /*OUTPUT how many steps did the optimization perform?               */
   int     t_limit;     /*INPUT maximum clocks to allow the process to run for
                         *      Note that this value will not be followed exaclty, there
                         *      will be overshoot                                            
                         *      Setting this value to zero allows infinite runtime          */
   int     step_limit;  /*INPUT maximum number of steps to allow the optimizer to perform   
                         *      Setting this value to zero allows infinite steps            */
   float   cutoff;      /*INPUT when the total change is below this value, exit. Default   
                         * value of 1e-4 will be used if this parameter is not set          */
   enum optimize_end_state 
           end_state;   /*OUTPUT returns why the optimization stopped                       */
} gsod;

/* produce nice printed output of the optimization data */
void print_graph_slam_optimization_data(gsod * data);

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
 * to exectue if the pose graph gets large. Note that
 * passing NULL to settings will result in defaults 
 * being used.                                           */
int optimize(pg * graph, gsod * settings);

#endif
