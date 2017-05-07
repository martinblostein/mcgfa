#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

// ETA UPDATES -----------

typedef struct {
    // struct used to pass extra arguments into optimization function
    int N;
    int p;
    double *mahalanobis;
    double *zbad;
} Info_eta;


double eta_objective_f(double eta, void *info_ptr) {
    // objective function that is maximized to produce eta update. From bottom of pg 25 of appendix
    int i;
    double output=0;

    Info_eta info = *((Info_eta*)info_ptr);
    for (i=0;i<info.N;i++) {
        output += info.zbad[i] * (info.p*log(eta) + info.mahalanobis[i]/eta);
    }
    return output;
}

void update_eta(double *eta, double eta_max, double *zbad, double **mahalanobis,
                    int N, int G, int p) {
    int i,g;
    double *zbad0 = malloc(sizeof(double)*N);

    Info_eta info;

    for(g=0; g<G; g++){
        for(i=0; i<N; i++) {
            zbad0[i] = zbad[i*G+g];
        }
        info.N = N; info.p = p; info.mahalanobis = mahalanobis[g]; info.zbad=zbad0;
        eta[g] = Brent_fmin(1, eta_max, &eta_objective_f, &info, 0.0001220703);
    }

    free(zbad0);
}



