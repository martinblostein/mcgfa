#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

// ---------- Z, V UPDATES ---------- //

// for CCC & CCU
void update_zv(double *log_dens, double *x, double *z, double *v, double *pi, double *max_log_dens,
               double log_c, double **mahal, double *eta, double *alpha, int N, int G, int p, int q) {

    int i,g;
    double e,d_alt=0;
    double *log_good_dens = malloc(sizeof(double)*N*G);
    double *log_bad_dens = malloc(sizeof(double)*N*G);
    double *log_dens_vec = malloc(sizeof(double)*G);

    for(i=0; i<N; i++){
        for(g=0; g<G; g++) {

            e = mahal[g][i]/2.0*(-1.0);

            log_good_dens[i*G+g] = e + log(pi[g]) - log_c;
            log_bad_dens[i*G+g] = e/eta[g] + log(pi[g]) - log_c - 0.5*p*log(eta[g]);
            log_dens[i*G+g] = log(alpha[g]*exp(log_good_dens[i*G+g]) + (1-alpha[g])*exp(log_bad_dens[i*G+g]));

        }

        if (G != 1) {
            for(g=0;g<G;g++) log_dens_vec[g] = log_dens[i*G+g];

            max_log_dens[i] = maximum_array(log_dens_vec,G);
            d_alt=0.0;

            for(g=0;g<G;g++) d_alt += exp(log_dens[i*G+g]-max_log_dens[i]);
        }


        for(g=0;g<G;g++) {
            if (G != 1) z[i*G+g] = exp(log_dens[i*G+g]-max_log_dens[i])/d_alt;
            v[i*G+g] = alpha[g]*exp(log_good_dens[i*G+g]-log_dens[i*G+g]);
            if (!isfinite(v[i*G+g]) || v[i*G+g]>1) v[i*G+g] = 1;
            if (v[i*G+g]<0) v[i*G+g] = 0;
        }
    }

    free(log_good_dens);
    free(log_bad_dens);
    free(log_dens_vec);

}

// for all other models
void update_zv2(double *log_dens, double *x, double *z, double *v, double *pi, double *max_log_dens,
                double *log_c, double **mahal, double *eta, double *alpha, int N, int G, int p, int q) {

    int i,g;
    double e,d_alt=0;
    double *log_good_dens = malloc(sizeof(double)*N*G);
    double *log_bad_dens = malloc(sizeof(double)*N*G);
    double *log_dens_vec = malloc(sizeof(double)*G);

    for(i=0; i<N; i++){
        for(g=0; g<G; g++) {

            e = mahal[g][i]/2.0*(-1.0);

            log_good_dens[i*G+g] = e + log(pi[g]) - log_c[g];
            log_bad_dens[i*G+g] = e/eta[g] + log(pi[g]) - log_c[g] - 0.5*p*log(eta[g]);
            log_dens[i*G+g] = log(alpha[g]*exp(log_good_dens[i*G+g]) + (1-alpha[g])*exp(log_bad_dens[i*G+g]));

        }

        if (G != 1) {
            for(g=0;g<G;g++) log_dens_vec[g] = log_dens[i*G+g];
            max_log_dens[i] = maximum_array(log_dens_vec,G);
            d_alt=0.0;
            for(g=0;g<G;g++) d_alt += exp(log_dens[i*G+g]-max_log_dens[i]);
        }

        for(g=0;g<G;g++) {
            if (G != 1) z[i*G+g] = exp(log_dens[i*G+g]-max_log_dens[i])/d_alt;
            v[i*G+g] = alpha[g]*exp(log_good_dens[i*G+g]-log_dens[i*G+g]);
            if (!isfinite(v[i*G+g]) || v[i*G+g]>1) v[i*G+g] = 1;
            if (v[i*G+g]<0) v[i*G+g] = 0;
        }
    }

    free(log_good_dens);
    free(log_bad_dens);
    free(log_dens_vec);

}

double maximum_array(double *array, int k) {
    // returns the maximum value of an array of doubles
    int i;
    double max = array[0];

    for(i=1;i<k;i++)
        if(array[i]> max) max=array[i];

        return max;
}
