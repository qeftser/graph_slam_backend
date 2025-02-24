
#ifndef SOLVER

#define SOLVER
#include "default.h"

/* some of the functions given here are sanity checks
 * and slower algorithms made as my understanding of
 * the situation evolved.
 */

/* prints the contents of a n*n square matrix */
void print_full_matrix(gsb_float * M, int n);

/* returns L, the lower triangular matrix for which L * Lt = A */
gsb_float * basic_cholesky_factorization(gsb_float * A, int n);

/* packed column representation of a matrix */
typedef struct packed_column_matrix {
   int n;       /* matrix order                      */
   int nnz;     /* number of non-zero elements       */
   int * colp;  /* array of start of column pointers */
   int * rx;    /* array of row indices              */
   gsb_float * val; /* array of off-digonal values       */
} pcmat;

typedef struct row_traversal_metadata {
   pcmat * A;  /* the matrix we are traversing               */
   int i;      /* the current row of the matrix              */
   int j;      /* the current column of the matrix           */
   int posij;  /* the current pos of the solver              */
   int nextj;   /* the next value of j                       */
   int * link; /* references to the next element in each row */
   int * pos;  /* position of the next element in each row   */
} rtmd;

typedef struct matrix_triplet {
   int row;   /* the row position of the value    */
   int col;   /* the column position of the value */
   gsb_float val; /* the actual value                 */
} mt;

/* auxillary method for comparing two matrix triplets */
int compare_matrix_triplets(const void * t1, const void * t2);
/* we input the list of values into a packed column format */
pcmat * construct_packed_column_matrix(mt * values, int n, int nnz);
/* get a uninitialized packed column matrix with values of the correct size */
pcmat * allocate_packed_column_matrix(int n, int nnz);
/* place sorted triplets into a pcmat that has already been initialized */
void load_packed_column_matrix(mt * values, pcmat * mat);
/* construct the data structure nessesary for traversing 
 * the packed matrix by row                              */
rtmd * get_row_traversal_metadata(pcmat * A);
/* perform a step on the row traversal. This will return a non-zero
 * value until a new value has been reached                         */
int row_traversal(rtmd * r);

/* construct the elimination tree for the pcmat A.
 * Also return the number of non-zero elements in the 
 * resulting L matrix. The elimination tree has the 
 * structure that each element in the column k has
 * the element parent[k] as it's parent, and each
 * element parent[k] has the element parent[parent[k]]
 * as it's parent. It can be shown that this linkage
 * eventually terminates at the value n                */
int construct_elimination_tree(pcmat * A);

/* sanity check */
void print_packed_column_matrix(pcmat * A);

/* forgot about these :/ */
void free_row_traversal_metadata(rtmd * r);
void free_packed_column_matrix(pcmat * A);

/* helper method that allows us to recursivly determine 
 * column structure                                    */
void merge_columns(pcmat * B, int j, int k, int A, int * ma);
/* helper method that transfers the merge array ma into
 * the column k of L
 * also resets the merge array                         */
void make_column(int k, int * ma, pcmat * L);
/* what we have been working up to! Bascially fill the
 * pcmat L with all the indicies that we will need to
 * perform the numerica factorization                 */
pcmat * compute_symbolic_factorization(pcmat * A);
/* self explanatory */
void perform_numerical_factorization(pcmat * A, pcmat * L);
/* solve the system of equations Ax = b using the factored
 * value L. We first solve for Ly = b, then for Ltx = y. 
 * The array b is overwritten with the values of x         */
void solve_system(pcmat * L, gsb_float * b);

#endif
