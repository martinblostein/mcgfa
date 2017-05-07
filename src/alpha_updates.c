#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

#define FMIN_TOL 0.0001220703
#define ALPHA_MAX 0.9995

typedef struct {
    int N;
    int p;
    double *z;
    double *v;
} Info_alpha;

double alpha_objective_f(double alpha, void *info_ptr) {
    int i;
    double output=0;

    Info_alpha info = *((Info_alpha*)info_ptr);
    for (i=0;i<info.N;i++) {
        output += info.z[i] * (info.v[i]*log(alpha) + (1-info.v[i])*log(1-alpha));
    }
    return -output;
}

void update_alpha_numerical(double *alpha, double alpha_min, double *z, double *v, int G, int N) {
    int i,g;

    double *z0 = malloc(sizeof(double)*N);
    double *v0 = malloc(sizeof(double)*N);

    Info_alpha info;

    for(g=0; g<G; g++){
        for(i=0; i<N; i++) {
            z0[i] = z[i*G+g];
            v0[i] = v[i*G+g];
        }

        info.N = N; info.z = z0; info.v=v0;

        alpha[g] = Brent_fmin(alpha_min,1, &alpha_objective_f, &info, FMIN_TOL);
    }
    free(z0); free(v0);
}

void update_alpha(double *alpha, double alpha_min, double *z, double *n, double *v, int G, int N) {
    int g, i;

    for (g=0; g<G; g++) {
        alpha[g] = 0;

        for (i=0; i<N; i++) {
            alpha[g] += (z[i*G+g] * v[i*G+g]);
        }
        alpha[g] = alpha[g]/n[g];

        if (alpha[g] < alpha_min) alpha[g] = alpha_min;
        if (alpha[g] > ALPHA_MAX) alpha[g] = ALPHA_MAX;
    }

}

void update_alpha_numerical_wrapper(double *alpha, double *alpha_min, double *z, double *v, int *G, int *N) {
    update_alpha_numerical(alpha, *alpha_min, z, v, *G, *N);
}

void update_alpha_wrapper(double *alpha, double *alpha_min, double *z, double *v, double *n, int *G, int *N) {
    update_alpha(alpha, *alpha_min, z, n, v, *G, *N);
}

