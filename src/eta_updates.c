#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

// ETA UPDATES -----------

void update_eta(double *eta, double eta_max, double *zbad, double **mahalanobis,
                int N, int G, int p, int fixed) {
    int i,g;
    double a, b, e;
    
    if (fixed) {
        a = 0;
        b = 0;
        for (g=0; g < G; g++) {
            for (i = 0; i < N; i++) {
                a += zbad[i*G+g];
                b += zbad[i*G+g] * mahalanobis[g][i];
            }
        }

        if (a == 0) {
            e = 1.01;
        } else {
            e = b / (p * a);
            
            if (e > eta_max) {
                e = eta_max;
            } else if (e < 1) {
                e = 1;
            }
        }

        for (g=0; g < G; g++) {
            eta[g] = e;
        }
        
        
    } else {
        
        for (g=0; g < G; g++) {
            a = 0;
            b = 0;
            for (i = 0; i < N; i++) {
                a += zbad[i*G+g];
                b += zbad[i*G+g] * mahalanobis[g][i];
            }
            
            if (a == 0) {
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
}
