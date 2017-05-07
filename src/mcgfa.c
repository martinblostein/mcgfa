#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"
#include "time.h"

double (*testfunc)(double*,double*,double*,int,int*,int,int,int,int,double*,
        double*,double*,double*,double*,double,double,double,int,int*);
typedef typeof(testfunc) funcType;
funcType funcs[8] = {aecm_CCC,aecm_CCU,aecm_CUC,aecm_CUU,aecm_UCC,aecm_UCU,aecm_UUC,aecm_UUU};

void mcgfa_c(double *x, double *z, double *v, double *bic, int *cls, int *q, int *p, int *G, int *N,
              double *mu, int *model, int *class_ind, double *lambda, double *psi, double *eta, double *alpha,
              double *eta_max, double *alpha_min, double *tol, int *max_it, int *iterations) {

    funcType aecm_func = funcs[*model-1];
    *bic = aecm_func(z, x, v, *class_ind, cls, *q, *p, *G, *N, mu, lambda, psi, eta, alpha, *eta_max, *alpha_min, *tol, *max_it, iterations);
}
