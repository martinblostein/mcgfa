#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<R.h>
#include "functions.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#endif


void generate_identity(int N, double *matrix){
    int i, j;
    for (i = 0; i < N; i ++)
        for (j = 0; j < N; j ++){
            matrix[i*N+j] = 0;
            if(i==j) matrix[i*N+i]=1;
        }
}

void GaussJordan(int N, double *MATRIX, double *INVERSE, double *det) {
    /* use generate_identity function to make INVERSE the identity. Then call
     * GaussJordan with MATRIX the matrix to invert. MATRIX gets turned into
     * the identity, so make a copy if you need to keep it. N is the dimension
     * of the matrix.*/
    int c, r, r_max, j, sign;
    double temp, v_v, v_max, factor = 0.0;
    generate_identity(N, INVERSE); det[0]=1; sign=0;

    /* Loop over all columns of A */
    for (c = 0; c < N; c++){
    /* Find row with the maximum value absolute value. */
        r_max = c;
        v_max = fabs(MATRIX[c*N+c]);
        for (r = c + 1; r < N; r++){
            v_v = fabs(MATRIX[c+r*N]);
            if (v_v > v_max){
                r_max = r;
                v_max = v_v;
            }
        }
 //First section
        /* Switch rows if necessary */
        if (r_max != c){
            for (j = c; j < N; j++){
                temp = MATRIX[c*N+j];
                MATRIX[c*N+j] = MATRIX[r_max*N+j];
                MATRIX[r_max*N+j] = temp;
            }
            for (j = 0; j < N; j++){
                temp = INVERSE[c*N+j];
                INVERSE[c*N+j] = INVERSE[r_max*N+j];
                INVERSE[r_max*N+j] = temp;
            }
            sign++;
        }
        /* Rescale current row so that diagonal element is 1 */
        factor = MATRIX[c*N+c];
        det[0]*= factor;
        /* printf("factor_Prod[%d]=%f\n",c,det[0]); */

        /*if (fabs(factor) == 0){
            printf("Matrix cannot be inverted.\n");
        }*/
        for (j = c; j < N; j++){
            MATRIX[c*N+j] /= factor;
        }
        for (j = 0; j < N; j++){
            INVERSE[c*N+j] /= factor;
        }
        /* Subtract current row from all rows below. */
        for (r = c + 1; r < N; r++){
            factor = MATRIX[c+r*N];
            for (j = c; j < N; j++){
                MATRIX[r*N+j] -= factor * MATRIX[c*N+j];
            }
            for (j = 0; j < N; j++)
                INVERSE[r*N+j] -= factor * INVERSE[c*N+j];
        }
    }

    if(sign%2 != 0) det[0] *= -1;

    /*loop over all rows from 1 to end*/
    for (c = 1; c < N; c++){
        /* Subtract current row from all rows above. */
        for (j = 0; j < N-c; j++){
            factor = MATRIX[N-c+j*N];
            for (r = 0; r < N; r++){
                MATRIX[j*N+r] =MATRIX[j*N+r] - factor * MATRIX[(N-c)*N+r];
                INVERSE[j*N+r] =INVERSE[j*N+r] - factor * INVERSE[(N-c)*N+r];
            }
        }
    }
}


void mx_mult(int m, int n, int q, double *a, double *b, double *r) {
    // Find matrix R = AB
    // Wrapper for FORTRAN function dgemm
    char notrans = 'N';
    double alpha = 1.0f;
    double beta = 0.0f;

    dgemm_(&notrans,&notrans, &q, &m, &n, &alpha, b, &q, a, &n, &beta, r, &q);
}


void mx_mult_diag1(int m, int n, double *a, double *b, double *r){
    // Find the diagonal of matrix R = AB
    int i, k;
    for(i=0; i<m; i++){
			r[i]=0.0;
            for(k=0; k<n; k++)
                r[i] += a[i*n+k]*b[i+k*m];
	}
}

/* Function to find r = bA */
void vec_mx_mult(int n, int q, double *a, double *b, double *r){
	int j,k;
        for(j=0; j<q; j++){
            r[j]=0.0;
            	for(k=0; k<n; k++){
                    r[j] += a[k]*b[j+k*q];
                }
        }
}

/* Function to find the transpose of an m*n matrix A; R=A' */
void mx_trans(int m, int n, double *a, double *r){
    int i,j;
    for(i=0; i<n; i++)
        for(j=0; j<m; j++)
            r[i*m+j]=a[i+j*n];
}
