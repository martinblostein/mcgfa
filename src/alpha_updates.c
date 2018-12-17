#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include "functions.h"

#define ALPHA_MAX 0.9995

void update_alpha(double *alpha, double alpha_min, double *z, double *n, double *v, int G, int N, int fixed) {
    int g, i;
    
    if (fixed) {
        double a = 0;
        
        for (g=0; g<G; g++) {
            for (i=0; i<N; i++) {
                a += (z[i*G+g] * v[i*G+g]);
            }
        }
        
        a = a/N;
        
        for (g=0; g<G; g++) {
            alpha[g] = a;
        }
    } else {
        for (g=0; g<G; g++) {
            alpha[g] = 0;
            
            for (i=0; i<N; i++) {
                alpha[g] += (z[i*G+g] * v[i*G+g]);
            }
            alpha[g] = alpha[g] / n[g];
            
            if (alpha[g] < alpha_min) alpha[g] = alpha_min;
            if (alpha[g] > ALPHA_MAX) alpha[g] = ALPHA_MAX;
        }
    }

}
