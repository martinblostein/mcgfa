#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

// ETA UPDATES -----------

void update_eta(double *eta, double eta_max, double *zbad, double **mahalanobis,
                           int N, int G, int p) {
    int i,g;
    double a, b;

    for (g=0; g < G; g++) {
        a = 0;
        b = 0;
        for (i = 0; i < N; i++) {
            a += zbad[i*G+g];
            b += zbad[i*G+g] * mahalanobis[g][i];
        }

        if (a == 0) {
            printf("Warning, all observations considered 'good' with probability 1 \
                        during eta update for group %d. Assigning eta = 1.1", g);
            eta[g] = 1.01;
        } else {
            eta[g] = b / (p * a);
        }

        if (eta[g] > eta_max) {
            eta[g] = eta_max;
        }  else if (eta[g] < 1) {
            eta[g] = 1;
        }
    }
}