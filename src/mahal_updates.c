#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

void update_mahal_CCC(double **mahal, double *x, double *lambda, double psi, double *mu, int N, int G, int p, int q) {
    int i, g, j;
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *sigma_inverse = malloc(sizeof(double)*p*p);

    // memory allocation for woodbury2 functions
    double *temp = malloc(sizeof(double)*q*p);
    double *temp2 = malloc(sizeof(double)*q*p);
    double *lvec = malloc(sizeof(double)*p);
    double *cvec = malloc(sizeof(double)*p);

    for(g=0; g<G; g++) {
        // calculate sigma inverse here
        for(j=0; j<p; j++) {
            mu0[j]= mu[g*p+j]; // g-th COLUMN of mu (mu_g)
        }

        woodburyA(lambda, psi, p, q, temp, temp2, sigma_inverse);

        for(i=0; i<N; i++) {
            for(j=0; j<p; j++) {
                x0[j] = x[i*p+j]; // i-th ROW of x

            }
            mahal[g][i] = woodburyB(x0, psi, mu0, p, q, lvec, cvec, sigma_inverse);
        }
    }

    free(x0); free(mu0); free(sigma_inverse);

    // memory deallocation for woodbury2 functions
    free(lvec); free(cvec); free(temp); free(temp2);

}

void update_mahal_CCU(double **mahal, double *x, double *lambda, double *psi, double *mu, int N, int G, int p, int q) {
    int i, g, j;
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *sigma_inverse = malloc(sizeof(double)*p*p);

    // memory allocation for woodbury2 functions
    double *temp = malloc(sizeof(double)*q*p);
    double *temp2 = malloc(sizeof(double)*q*p);
    double *lvec = malloc(sizeof(double)*p);
    double *cvec = malloc(sizeof(double)*p);

    for(g=0; g<G; g++) {
        // calculate sigma inverse here
        for(j=0; j<p; j++) {
            mu0[j]= mu[g*p+j]; // g-th COLUMN of mu (mu_g)
        }

        woodbury2A(lambda, psi, p, q, temp, temp2, sigma_inverse);

        for(i=0; i<N; i++) {
            for(j=0; j<p; j++) {
                x0[j] = x[i*p+j]; // i-th ROW of x

            }
            mahal[g][i] = woodbury2B(x0, psi, mu0, p, q, lvec, cvec, sigma_inverse);
        }
    }

    free(x0); free(mu0); free(sigma_inverse);

    // memory deallocation for woodbury2 functions
    free(lvec); free(cvec); free(temp); free(temp2);
}

void update_mahal_CUC(double **mahal, double *x, double *lambda, double *psi, double *mu, int N, int G, int p, int q) {
    int i, g, j;
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *sigma_inverse = malloc(sizeof(double)*p*p);

    // memory allocation for woodbury2 functions
    double *temp = malloc(sizeof(double)*q*p);
    double *temp2 = malloc(sizeof(double)*q*p);
    double *lvec = malloc(sizeof(double)*p);
    double *cvec = malloc(sizeof(double)*p);

    for(g=0; g<G; g++) {
        // calculate sigma inverse here
        for(j=0; j<p; j++) {
            mu0[j]= mu[g*p+j]; // g-th COLUMN of mu (mu_g)
        }

        woodburyA(lambda, psi[g], p, q, temp, temp2, sigma_inverse);

        for(i=0; i<N; i++) {
            for(j=0; j<p; j++) {
                x0[j] = x[i*p+j]; // i-th ROW of x

            }
            mahal[g][i] = woodburyB(x0, psi[g], mu0, p, q, lvec, cvec, sigma_inverse);
        }
    }

    free(x0); free(mu0); free(sigma_inverse);

    // memory deallocation for woodbury2 functions
    free(lvec); free(cvec); free(temp); free(temp2);
}

void update_mahal_CUU(double **mahal, double *x, double *lambda, double *psi, double *mu, int N, int G, int p, int q) {
    int i, g, j;
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *psi0 = malloc(sizeof(double)*p);
    double *sigma_inverse = malloc(sizeof(double)*p*p);

    // memory allocation for woodbury2 functions
    double *temp = malloc(sizeof(double)*q*p);
    double *temp2 = malloc(sizeof(double)*q*p);
    double *lvec = malloc(sizeof(double)*p);
    double *cvec = malloc(sizeof(double)*p);

    for(g=0; g<G; g++) {
        // calculate sigma inverse here
        for(j=0; j<p; j++) {
            psi0[j]= psi[g*p+j]; // psi vector for group g
            mu0[j]= mu[g*p+j]; // g-th COLUMN of mu (mu_g)
        }

        woodbury2A(lambda, psi0, p, q, temp, temp2, sigma_inverse);

        for(i=0; i<N; i++) {
            for(j=0; j<p; j++) {
                x0[j] = x[i*p+j]; // i-th ROW of x

            }
            mahal[g][i] = woodbury2B(x0, psi0, mu0, p, q, lvec, cvec, sigma_inverse);
        }
    }

    free(x0); free(mu0); free(psi0); free(sigma_inverse);

    // memory deallocation for woodbury2 functions
    free(lvec); free(cvec); free(temp); free(temp2);

}

void update_mahal_UCC(double **mahal, double *x, double **lambda, double psi, double *mu, int N, int G, int p, int q) {
    int i, g, j;
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *sigma_inverse = malloc(sizeof(double)*p*p);

    // memory allocation for woodbury2 functions
    double *temp = malloc(sizeof(double)*q*p);
    double *temp2 = malloc(sizeof(double)*q*p);
    double *lvec = malloc(sizeof(double)*p);
    double *cvec = malloc(sizeof(double)*p);

    for(g=0; g<G; g++) {
        // calculate sigma inverse here
        for(j=0; j<p; j++) {
            mu0[j]= mu[g*p+j]; // g-th COLUMN of mu (mu_g)
        }

        woodburyA(lambda[g], psi, p, q, temp, temp2, sigma_inverse);

        for(i=0; i<N; i++) {
            for(j=0; j<p; j++) {
                x0[j] = x[i*p+j]; // i-th ROW of x

            }
            mahal[g][i] = woodburyB(x0, psi, mu0, p, q, lvec, cvec, sigma_inverse);
        }
    }

    free(x0); free(mu0); free(sigma_inverse);

    // memory deallocation for woodbury2 functions
    free(lvec); free(cvec); free(temp); free(temp2);
}

void update_mahal_UCU(double **mahal, double *x, double **lambda, double *psi, double *mu, int N, int G, int p, int q) {
    int i, g, j;
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *sigma_inverse = malloc(sizeof(double)*p*p);

    // memory allocation for woodbury2 functions
    double *temp = malloc(sizeof(double)*q*p);
    double *temp2 = malloc(sizeof(double)*q*p);
    double *lvec = malloc(sizeof(double)*p);
    double *cvec = malloc(sizeof(double)*p);

    for(g=0; g<G; g++) {
        // calculate sigma inverse here
        for(j=0; j<p; j++) {
            mu0[j]= mu[g*p+j]; // g-th COLUMN of mu (mu_g)
        }

        woodbury2A(lambda[g], psi, p, q, temp, temp2, sigma_inverse);

        for(i=0; i<N; i++) {
            for(j=0; j<p; j++) {
                x0[j] = x[i*p+j]; // i-th ROW of x

            }
            mahal[g][i] = woodbury2B(x0, psi, mu0, p, q, lvec, cvec, sigma_inverse);
        }
    }

    free(x0); free(mu0);  free(sigma_inverse);

    // memory deallocation for woodbury2 functions
    free(lvec); free(cvec); free(temp); free(temp2);
}

void update_mahal_UUC(double **mahal, double *x, double **lambda, double *psi, double *mu, int N, int G, int p, int q) {
    int i, g, j;
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *sigma_inverse = malloc(sizeof(double)*p*p);

    // memory allocation for woodbury2 functions
    double *temp = malloc(sizeof(double)*q*p);
    double *temp2 = malloc(sizeof(double)*q*p);
    double *lvec = malloc(sizeof(double)*p);
    double *cvec = malloc(sizeof(double)*p);

    for(g=0; g<G; g++) {
        // calculate sigma inverse here
        for(j=0; j<p; j++) {
            mu0[j]= mu[g*p+j]; // g-th COLUMN of mu (mu_g)
        }

        woodburyA(lambda[g], psi[g], p, q, temp, temp2, sigma_inverse);

        for(i=0; i<N; i++) {
            for(j=0; j<p; j++) {
                x0[j] = x[i*p+j]; // i-th ROW of x

            }
            mahal[g][i] = woodburyB(x0, psi[g], mu0, p, q, lvec, cvec, sigma_inverse);
        }
    }

    free(x0); free(mu0); free(sigma_inverse);

    // memory deallocation for woodbury2 functions
    free(lvec); free(cvec); free(temp); free(temp2);
}

void update_mahal_UUU(double **mahal, double *x, double **lambda, double *psi, double *mu, int N, int G, int p, int q) {
    int i, g, j;
    double *x0 = malloc(sizeof(double)*p);
    double *mu0 = malloc(sizeof(double)*p);
    double *psi0 = malloc(sizeof(double)*p);
    double *sigma_inverse = malloc(sizeof(double)*p*p);

    // memory allocation for woodbury2 functions
    double *temp = malloc(sizeof(double)*q*p);
    double *temp2 = malloc(sizeof(double)*q*p);
    double *lvec = malloc(sizeof(double)*p);
    double *cvec = malloc(sizeof(double)*p);


    for(g=0; g<G; g++) {
        // calculate sigma inverse here
        for(j=0; j<p; j++) {
            psi0[j]= psi[g*p+j]; // psi vector for group g
            mu0[j]= mu[g*p+j]; // g-th COLUMN of mu (mu_g)
        }

        woodbury2A(lambda[g], psi0, p, q, temp, temp2, sigma_inverse);

        for(i=0; i<N; i++) {
            for(j=0; j<p; j++) {
                x0[j] = x[i*p+j]; // i-th ROW of x

            }
            mahal[g][i] = woodbury2B(x0, psi0, mu0, p, q, lvec, cvec, sigma_inverse);
        }
    }

    free(x0); free(mu0); free(psi0); free(sigma_inverse);

    // memory deallocation for woodbury2 functions
    free(lvec); free(cvec); free(temp); free(temp2);
}

// --------------------------------------------------------------------------------------

void woodburyA(double *lambda, double psi, int p, int q, double *A, double *B, double *C) {
    // This function computes:
    // SIGMA^{-1} = LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    // and stores it in the array C. The arrays A & B are temporary containers only.
    // --> same as woodburyA but takes a vector argument for psi

    int i,j;
    double det;

    // A <- LAMBDA'
    mx_trans(p, q, lambda, A);
    // A <- LAMBDA' * PSI^{-1}
    for(i=0;i<q;i++) for(j=0;j<p;j++) A[i*p+j] /= psi;
    // B <- LAMBDA' * PSI^{-1} * LAMBDA
    mx_mult(q, p, q, A, lambda, B);
    // B <- I + LAMBDA' * PSI^{-1} * LAMBDA
    for(i=0;i<q;i++) B[i*q+i] += 1.0;
    // C <- (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1}
    GaussJordan(q, B, C, &det);
    // A <- LAMBDA'
    mx_trans(p, q, lambda, A);
    // B <- (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    mx_mult(q, q, p, C, A, B);
    // C <- LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    mx_mult(p, q, p, lambda, B, C);
    // C <- PSI^{-1} * LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    for(i=0;i<p;i++) for(j=0;j<p;j++) C[i*p+j] /= psi;
    // C <- PSI^{-1} * LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA' * PSI^{-1}
    for(i=0;i<p;i++) for(j=0;j<p;j++) C[i*p+j] /= psi;

}

double woodburyB(double *x, double psi, double *mu, int p, int q,
                 double *lvec, double *cvec, double *cp) {
    // This function takes SIGMA^{-1} as calculated by woodbury2A and uses it to compute
    // the Mahalanobis distance of x, given SIGMA{-1} and psi.
    // --> same as woodburyB but takes a vector argument for psi

    int j;
    double lhs=0, rhs=0;

    // lhs <- (x-mu)' * PSI^{-1} * (x-mu)
    for (j=0;j<p;j++) lhs += (x[j]-mu[j])*(x[j]-mu[j])/psi;

    // lvec <- (x-mu)'
    for (j=0;j<p;j++) lvec[j] = (x[j]-mu[j]);
    // cvec <- (x-mu)' * PSI^{-1} * LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    vec_mx_mult(p,p,lvec,cp,cvec);
    // rhs <- (x-mu)' * PSI^{-1} * LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA' * PSI^{-1} * (x-mu)
    for (j=0;j<p;j++) rhs += cvec[j]*(x[j]-mu[j]);

    // (lsh - rhs) is thus the mahalanobis distance of x_i as a member of the gth group
    return (lhs-rhs);
}

void woodbury2A(double *lambda, double *psi, int p, int q, double *A, double *B, double *C) {
    // This function computes:
    // SIGMA^{-1} = LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    // and stores it in the array C. The arrays A & B are temporary containers only.
    // --> same as woodburyA but takes a vector argument for psi

    int i,j;
    double det;

    // A <- LAMBDA'
    mx_trans(p, q, lambda, A);
    // A <- LAMBDA' * PSI^{-1}
    for(i=0;i<q;i++) for(j=0;j<p;j++) A[i*p+j] /= psi[j];
    // B <- LAMBDA' * PSI^{-1} * LAMBDA
    mx_mult(q, p, q, A, lambda, B);
    // B <- I + LAMBDA' * PSI^{-1} * LAMBDA
    for(i=0;i<q;i++) B[i*q+i] += 1.0;
    // C <- (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1}
    GaussJordan(q, B, C, &det);
    // A <- LAMBDA'
    mx_trans(p, q, lambda, A);
    // B <- (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    mx_mult(q, q, p, C, A, B);
    // C <- LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    mx_mult(p, q, p, lambda, B, C);
    // C <- PSI^{-1} * LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    for(i=0;i<p;i++) for(j=0;j<p;j++) C[i*p+j] /= psi[i];
    // C <- PSI^{-1} * LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA' * PSI^{-1}
    for(i=0;i<p;i++) for(j=0;j<p;j++) C[i*p+j] /= psi[j];

}

double woodbury2B(double *x, double *psi, double *mu, int p, int q,
                  double *lvec, double *cvec, double *cp) {
    // This function takes SIGMA^{-1} as calculated by woodbury2A and uses it to compute
    // the Mahalanobis distance of x, given SIGMA{-1} and psi.
    // --> same as woodburyB but takes a vector argument for psi

    int j;
    double lhs=0, rhs=0;

    // lhs <- (x-mu)' * PSI^{-1} * (x-mu)
    for (j=0;j<p;j++) lhs += (x[j]-mu[j])*(x[j]-mu[j])/psi[j];

    // lvec <- (x-mu)'
    for (j=0;j<p;j++) lvec[j] = (x[j]-mu[j]);
    // cvec <- (x-mu)' * PSI^{-1} * LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA'
    vec_mx_mult(p,p,lvec,cp,cvec);
    // rhs <- (x-mu)' * PSI^{-1} * LAMBDA * (I + LAMBDA' * PSI^{-1} * LAMBDA)^{-1} * LAMBDA' * PSI^{-1} * (x-mu)
    for (j=0;j<p;j++) rhs += cvec[j]*(x[j]-mu[j]);

    // (lsh - rhs) is thus the mahalanobis distance of x_i as a member of the gth group
    return (lhs-rhs);
}
