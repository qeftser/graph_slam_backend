
#include "solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void print_full_matrix(float * M, int n) {
   for (int i = 0; i < n; ++i) {
      printf("[ ");
      for (int j = 0; j < n; ++j) {
         printf("%+03.3f ",M[j+(i*n)]);
      }
      printf("]\n");
   }
}

float * basic_cholesky_factorization(float * A, int n) {
   /* L is our return matrix */
   float * L = calloc(sizeof(float),n*n);

   /* copy the lower triangular section of A into L */
   for (int i = 0; i < n; ++i) {
      for (int j = 0; j <= i; ++j) {
         L[j+(i*n)] = A[j+(i*n)];
      }
   }

   /* solve the factorization */
   for (int k = 0; k < n; ++k) {
      for (int j = 0; j < k; ++j) {
         for (int i = k; i < n; ++i) {
            L[k+(i*n)] -= L[j+(k*n)]*L[j+(i*n)];
         }
      }
      L[k+(k*n)] = sqrt(L[k+(k*n)]);
      for (int i = k+1; i < n; ++i) {
         L[k+(i*n)] /= L[k+(k*n)];
      }
   }

   return L;
}

/* assume the user is not dumb enough to add the same element twice :/ */
/* I am the user so this makes me nervous                              */
int compare_matrix_triplets(const void * t1, const void * t2) {
   if (((mt *)t1)->col < ((mt *)t2)->col)
      return -1;
   else if (((mt *)t1)->col > ((mt *)t2)->col)
      return 1;
   else {
      if (((mt *)t1)->row < ((mt *)t2)->row)
         return -1;
      return 1;
   }
   return 0;
}

/* Note: this function assumes that there will be at least
 * one non-zero value at each column in the matrix         */
pcmat * construct_packed_column_matrix(mt * values, int n, int nnz) {
   /* allocate memory */
   pcmat * ret = calloc(sizeof(pcmat),1);
   ret->colp = calloc(sizeof(int),n+1);
   ret->rx = calloc(sizeof(int),nnz);
   ret->val = calloc(sizeof(float),nnz);

   /* fill basic values */
   ret->n = n;
   ret->nnz = nnz;

   /* insert triples into the structure */
   int col = 0, row = 0, last_col = -1;
   qsort(values,nnz,sizeof(mt),compare_matrix_triplets);
   for (int i = 0; i < nnz; ++i) {
      ret->rx[row] = values[i].row;
      ret->val[row] = values[i].val;
      if (values[i].col != last_col) {
         last_col = values[i].col;
         ret->colp[col] = row;
         ++col;
      }
      ++row; 
   }
   ret->colp[col] = row;

   return ret;
}

pcmat * allocate_packed_column_matrix(int n, int nnz) {
   pcmat * ret = calloc(sizeof(pcmat),1);
   ret->colp = calloc(sizeof(int),n+1);
   ret->rx = calloc(sizeof(int),nnz);
   ret->val = calloc(sizeof(float),nnz);
   ret->n = n;
   ret->nnz = nnz;
   return ret;
}

/* same as the second half of construct_packed_column_matrix */
void load_packed_column_matrix(mt * values, pcmat * mat) {
   int col = 0, row = 0, last_col = -1;
   qsort(values,mat->nnz,sizeof(mt),compare_matrix_triplets);
   for (int i = 0; i < mat->nnz; ++i) {
      mat->rx[row] = values[i].row;
      mat->val[row] = values[i].val;
      if (values[i].col != last_col) {
         last_col = values[i].col;
         mat->colp[col] = row;
         ++col;
      }
      ++row;
   }
   mat->colp[col] = row;
}

void print_packed_column_matrix(pcmat * A) {
   printf("n: %5d nnz: %5d\n",A->n,A->nnz);
   int col = 0;
   int pos = 0;
   /* print columns. should be a do-while xD */
   printf("[ ");
step:
   if (pos == A->colp[col]) {
      printf("%6d ",A->colp[col]);
      ++col;
   }
   else
      printf("       ");
   ++pos;
   if (pos < A->nnz)
      goto step;
   printf("]\n");

   /* print rows */
   printf("[ ");
   for (int i = 0; i < A->nnz; ++i)
      printf("%6d ",A->rx[i]);
   printf("]\n");

   /* print values */
   printf("[ ");
   for (int i = 0; i < A->nnz; ++i)
      printf("%+06.2f ",A->val[i]);
   printf("]\n");
}

rtmd * get_row_traversal_metadata(pcmat * A) {
   /* allocate memory */
   rtmd * ret = calloc(sizeof(rtmd),1);
   ret->link = malloc(sizeof(int)*A->n);
   ret->pos = malloc(sizeof(int)*A->n);
   for (int i = 0; i < A->n; ++i)
      ret->link[i] = ret->pos[i] = -1;

   /* set initial values */
   ret->i = -1;
   ret->j = -1;
   ret->posij = -1;
   ret->nextj = -1;
   ret->A = A;

   /* a note for why this deviates from the paper - 
    * because we have to construct the traversal value
    * as opposed to using static variables we can do away
    * with using i = -1 as a flag. Instead we do our
    * initialization here. Also this stuff actually ends
    * up being way more elegant if we 1 index our arrays :/ */

   return ret;
}

/* mad respect to the guy that produced this */
int row_traversal(rtmd * r) {
   if (r->j == -1) {
      r->i = r->i + 1;
      r->j = r->i;
      r->posij = r->A->colp[r->i];
   }
   else {
      r->j = r->nextj;
      if (r->j == -1) return r->j;
      r->posij = r->pos[r->j];
   }
   r->nextj = r->link[r->j];
   r->link[r->j] = -1;
   int nextdown = r->posij + 1;
   if (nextdown < r->A->colp[r->j + 1]) {
      r->pos[r->j] = nextdown;
      int id = r->A->rx[nextdown];
      r->link[r->j] = r->link[id];
      r->link[id] = r->j;
   }
   return r->j;
}

void free_row_traversal_metadata(rtmd * r) {
   free(r->link); free(r->pos); free(r);
}

void free_packed_column_matrix(pcmat * A) {
   free(A->colp); free(A->rx); free(A->val); free(A);
}

int construct_elimination_tree(pcmat * A) {
   /* initialize values */
   int nnz = 0;
   int * touched = malloc(sizeof(int)*A->n);
   int * parent = malloc(sizeof(int)*A->n);
   for (int i = 0; i < A->n; ++i)
      touched[i] = parent[i] = -1;

   /* traverse rows of A */
   rtmd * r = get_row_traversal_metadata(A);
   for (int i = 0; i < A->n; ++i) {
      while (row_traversal(r) != -1) {
         if (r->i == r->j) {
            ++nnz;
            touched[r->j] = r->i;
         }
         else {
            int js = r->j;
            while (touched[js] != r->i) {
               touched[js] = r->i;
               ++nnz;
               if (parent[js] == -1) {
                  parent[js] = r->i;
                  break;
               }
               js = parent[js];
            }
         }
      }
   }
   free_row_traversal_metadata(r);
   free(parent);
   free(touched);

   return nnz;
}

/* basically just an extremly dense set union with 
 * some implicit stacks                            */
void merge_columns(pcmat * B, int j, int k, int A, int * ma) {
   int m = k;
   /* loop over elements in column j of B */
   for (int ii = B->colp[j] + 1+A; ii < B->colp[j+1]; ++ii) {
      int i = B->rx[ii];
      /* search for m and m1 with m < i <= m1; */
      int m1 = m;
      while (i > m1) {
         m = m1;
         m1 = ma[m];
      }
      if (i != m1) {
         /* insert i in ma */
         ma[m] = i;
         ma[i] = m1;
      }
      m = i;
   }
}

void make_column(int k, int * ma, pcmat * L) {
   if (k == 0) 
      L->colp[0] = 0; 
   int ii = L->colp[k];
   int m = k;
   while (m < L->n) {
      L->rx[ii] = m;
      ++ii;
      int mt = ma[m];
      ma[m] = L->n;
      m = mt;
   }
   L->colp[k+1] = ii;
}

pcmat * compute_symbolic_factorization(pcmat * A) {
   /* initialize L */
   pcmat * L = malloc(sizeof(pcmat));
   int nnz = construct_elimination_tree(A);
   L->nnz = nnz;
   L->n = A->n;
   L->colp = malloc(sizeof(int)*(L->n+1));
   L->rx = malloc(sizeof(int)*nnz);
   L->val = malloc(sizeof(float)*nnz);

   /* initialize babysitter and merge arrays */
   int * bs = malloc(sizeof(int)*L->n);
   int * ma = malloc(sizeof(int)*L->n);
   for (int i = 0; i < L->n; ++i) {
      bs[i] = -1;
      ma[i] = L->n;
   }

   /* main loop on columns of A */
   int j, jt;
   for (int k = 0; k < A->n; ++k) {
      merge_columns(A,k,k,0,ma);
      j = bs[k];
      while (j != -1) {
         merge_columns(L,j,k,1,ma);
         jt = bs[j]; bs[j] = -1; j = jt;
      }
      /* set up the kth column of L */
      make_column(k,ma,L);
      /* update the babysitter */
      if (k != L->n-1) {
         j = L->rx[L->colp[k] + 1]; /* j is the parent of k */
         while (j != -1) {
            jt = j; j = bs[j];
         }
         bs[jt] = k;
      }
   }
   free(ma);
   free(bs);

   return L;
}

void perform_numerical_factorization(pcmat * A, pcmat * L) {

   /* an extra small value is needed to offset issues
    * with numerical stability */
   static float epsilion = 0.0000001;

   /* setup */
   rtmd * r = get_row_traversal_metadata(L);
   float * accum = calloc(sizeof(float),A->n);
   int ii, kx, i;
   float Lkj, Lkkinv;


   for (kx = 0; kx < L->n; ++kx) { /* process column k */
      while (row_traversal(r) != -1) {
         if (r->i == r->j) { /* initialize accum */
            for (ii = L->colp[r->i]; ii < L->colp[r->i+1]; ++ii)
               accum[L->rx[ii]] = 0.0;
            for (ii = A->colp[r->i]; ii < A->colp[r->i+1]; ++ii)
               accum[A->rx[ii]] = A->val[ii];
         }
         else { /* subtract L[k:n,j] from L[k:n,k] */
            Lkj = L->val[r->posij];
            for (ii = r->posij; ii < L->colp[r->j+1]; ++ii) {
               i = L->rx[ii];
               accum[i] = accum[i] - Lkj*L->val[ii];
            }
         }
      }
      /* move L[k:n,k] from accum to L, adjusting it's components */
      for (ii = L->colp[r->i]; ii < L->colp[r->i+1]; ++ii) {
         i = L->rx[ii];
         if (i == r->i) {
            L->val[ii] = sqrt(accum[i] + epsilion);
            Lkkinv = 1.0 / L->val[ii];
         }
         else {
            L->val[ii] = Lkkinv * accum[i];
         }
      }
   }
   free(accum);
   free_row_traversal_metadata(r);
}

void solve_system(pcmat * L, float * b) {

   /* solve for Ly = b */
   for (int j = 0; j < L->n; ++j) {
      b[j] = b[j] / L->val[L->colp[j]];
      for (int ii = L->colp[j] + 1; ii < L->colp[j+1]; ++ii) {
         int i = L->rx[ii];
         b[i] = b[i] - b[j] * L->val[ii];
      }
   }

   /* solve for Ltx = y */
   for (int j = L->n - 1; j >= 0; --j) {
      for (int ii = L->colp[j]+1; ii < L->colp[j+1]; ++ii) {
         int i = L->rx[ii];
         b[j] = b[j] - b[i] * L->val[ii];
      }
      b[j] = b[j] / L->val[L->colp[j]];
   }

   /* the solution x is now contained in b */
}
