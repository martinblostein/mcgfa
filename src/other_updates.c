#include <R.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "functions.h"

// ---------- PROPORTION UPDATES ----------//

void update_n(double *n, double *z, int G, int N){
    int g,i;
    for(g=0; g<G; g++){
        n[g]=0.0;
        for(i=0; i<N; i++){
            n[g] += z[g+i*G];
        }
    }
}

void update_pi(double *pi, double *n, int G, int N){
    int g;
    for(g=0; g<G; g++) pi[g] = n[g]/N;
}

void update_correction(double* correction, double* v, double* eta, int G, int N) {
    // new
    int i, g;

    for (g=0;g<G;g++) {
        for (i=0;i<N;i++) {
            correction[i*G+g] = v[i*G+g] + (1 - v[i*G+g])/eta[g];
        }
    }

}

// ---------- MU UPDATE ---------- //

void update_mu(double *mu, double *n, double *x, double *z, double *correction, int G, int N, int p) {
    // mu is
    // modified to incorporate correction
    int i,j,g;

    double *numerator = malloc(sizeof(double)*G*p);
    double *denominator = malloc(sizeof(double)*G*p);

    for(g=0; g<G; g++){

        for(j=0; j<p; j++){
            numerator[g*p+j]=0.0;
            denominator[g*p+j]=0.0;
            for(i=0; i<N; i++) {
                numerator[g*p+j] += z[i*G+g] * correction[i*G+g] * x[i*p+j];
                denominator[g*p+j] += z[i*G+g] * correction[i*G+g];
            }
            mu[g*p+j] = numerator[g*p+j]/denominator[g*p+j];
        }
    }
    
    free(numerator);
    free(denominator);
}

// ---------- S UPDATES ----------//

void update_stilde(double *sampcovtilde, double *x, double *z, double *mu, double *correction,
                   int G, int N, int p){
    int i,j,k,g;

    for(j=0; j<p; j++){
        for(k=0; k<p; k++){
            sampcovtilde[j*p+k]=0.0;
            for(g=0;g<G;g++)
                for(i=0; i<N; i++){

                    sampcovtilde[j*p+k] += z[g+i*G] * correction[g+i*G] * (x[j+i*p] - mu[j+g*p]) * (x[k+i*p]-mu[k+g*p]);
                }
                sampcovtilde[j*p+k]/=N;
        }
    }
}

void update_sg(double **S, double *x, double *z, double *mu, double *correction, double *n, int p, int G, int N){

    int i,j,k,g;

    for(g=0;g<G;g++){

        for(j=0; j<p; j++){
            for(k=j; k<p; k++) {
                S[g][j*p+k]=0.0;
                for(i=0; i<N; i++){
                    S[g][j*p+k] += z[g+i*G] * correction[g+i*G] * (x[j+i*p]-mu[g*p+j])*(x[k+i*p]-mu[g*p+k]);
                }
                S[g][j*p+k] /= n[g];
            }
            for(k=0; k<j; k++) {
                S[g][j*p+k] = S[g][k*p+j];
            }
        }
    }
}

// ---------- BETA UPDATES ---------- //

void update_beta_iso(double *beta, double psi, double *lambda, int p, int q){
    int i,j;
    double det[1];

    double *lhs = malloc(sizeof(double)*q*p);
    double *rhs = malloc(sizeof(double)*p*p);
    double *cp = malloc(sizeof(double)*q*q);
    double *result = malloc(sizeof(double)*p*p);


    /* Lambda'/psi */
    mx_trans(p, q, lambda, lhs);
    for(i=0;i<q;i++)
        for(j=0;j<p;j++)
            lhs[i*p+j]/=psi; /*psi[j]*/
    /* Lambda'/psi* Lambda */
    mx_mult(q, p, q, lhs, lambda, cp);

    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++){
        for(j=0;j<q;j++){
            result[i*q+j] = cp[i*q+j];
            if (i==j) result[i*q+i] += 1.0;
        }
    }
    GaussJordan(q, result, rhs, det);
    /* work out rhs*/
    mx_mult(q,q,q,cp,rhs,result);
    mx_mult(q,q,p,result,lhs,rhs);

    for(i=0;i<q;i++){
        for(j=0;j<p;j++){
            beta[i*p+j] = lhs[i*p+j] - rhs[i*p+j];
        }
    }
    free(lhs); free(result); free(rhs); free(cp);

}

void update_beta_niso(double *beta, double *Psi, double *lambda, int p, int q){
    int i,j;
    double det[1];
    double *lhs = malloc(sizeof(double)*p*p);
    double *rhs = malloc(sizeof(double)*p*p);
    double *cp = malloc(sizeof(double)*p*p);
    double *result = malloc(sizeof(double)*p*p);

    /* Lambda'/psi */
    mx_trans(p, q, lambda, lhs);
    for(i=0;i<q;i++)
        for(j=0;j<p;j++)
            lhs[i*p+j]/=Psi[j];
    /* Lambda'/psi* Lambda */
    mx_mult(q, p, q, lhs, lambda, cp);

    /* (I + Lambda'/psi*Lambda)^{-1} */
    for(i=0; i<q;i++){
        for(j=0;j<q;j++){
            result[i*q+j] = cp[i*q+j];
            if (i==j) result[i*q+i] += 1.0;
        }
    }
    GaussJordan(q, result, rhs, det);
    /* work out rhs*/
    mx_mult(q,q,q,cp,rhs,result);
    mx_mult(q,q,p,result,lhs,rhs);

    for(i=0;i<q;i++){
        for(j=0;j<p;j++){
            beta[i*p+j] = lhs[i*p+j] - rhs[i*p+j];
        }
    }
    free(lhs); free(result); free(rhs); free(cp);

}

// ---------- THETA UPDATE ---------- //

void update_theta(double *theta, double *beta, double *lambda, double *sampcovtilde, int p, int q){
    int i,j;

    double *tmp = malloc(sizeof(double)*p*p);
    double *r_1 = malloc(sizeof(double)*q*q);
    double *r_2 = malloc(sizeof(double)*q*p);
    double *r_3 = malloc(sizeof(double)*q*q);
    double *id_q = malloc(sizeof(double)*q*q);
    generate_identity(q, id_q);

    /* Work out beta*lambda */
    mx_mult(q, p, q, beta, lambda, r_1);

    /* Now, work out beta*sampcovtilde*beta' */
    mx_mult(q, p, p, beta, sampcovtilde, r_2);
    mx_trans(q, p, beta, tmp);
    mx_mult(q, p, q, r_2, tmp, r_3);

    for(i=0;i<q;i++){
        for(j=0;j<q;j++){
            theta[i*q+j] = id_q[i*q+j]-r_1[i*q+j]+r_3[i*q+j];
        }
    }
    free(id_q); free(tmp); free(r_1); free(r_2); free(r_3);

}

// The below functions are essentially unchanged from pgmm:

// ---------- LAMBDA UPDATES ---------- //

void update_lambda(double *lambda, double *beta, double *s, double *theta, int p, int q){
    // update lambda for all models except CCU & CUC
    int i,j;
    double  det[1];

    double *tran = malloc(sizeof(double)*p*q);
    double *res1 = malloc(sizeof(double)*p*q);
    double *res2 = malloc(sizeof(double)*q*q);
    double *res3 = malloc(sizeof(double)*q*q);

    // sampcovtilde * beta'
    mx_trans(q,p,beta,tran);
    mx_mult(p,p,q,s,tran,res1);

    // Copy of theta ahead of Gauss-Jordan
    // (Second argument of GJ function becomes identity.)
    for(i=0;i<q;i++)
        for(j=0;j<q;j++)
            res3[i*q+j] = theta[i*q+j];

    /* lambda=sampcovtilde*beta'*theta^{-1} */
    GaussJordan(q, res3, res2, det);
    mx_mult(p,q,q,res1,res2,lambda);
    free(tran); free(res1); free(res2); free(res3);

}

void update_lambda_cuu(double *lambda, double **beta, double **s,
                       double **theta, double *n, double *Psi, int p, int q, int G){

    int i,j,k,g;
    double det[1];
    double *tran = malloc(sizeof(double)*p*q);
    double *res1 = malloc(sizeof(double)*p*q);
    double *res2 = malloc(sizeof(double)*p*q);
    double *res3 = malloc(sizeof(double)*q*q);
    double *result = malloc(sizeof(double)*q);
    double *lambda0 = malloc(sizeof(double)*q);


    /* Compute the RHS --- only needs to happen once */
    // sum_g { sampcov_g * beta_g' * n_g/psi_g }
    for(g=0;g<G;g++){
        mx_trans(q,p,beta[g],tran);
        mx_mult(p,p,q,s[g],tran,res1);
        if (g==0) {
            for(i=0;i<p;i++){
                for(j=0;j<q;j++){
                    res2[i*q+j] = res1[i*q+j]*n[g]/Psi[g*p+i];
                }
            }
        } else {
            for(i=0;i<p;i++)
                for(j=0;j<q;j++)
                    res2[i*q+j] += res1[i*q+j]*n[g]/Psi[g*p+i];
        }
    }


    /* Now solve for lambda row-by-row */
    for (i=0;i<p;i++){
        /* First, compute the theta sum*/
        for(g=0;g<G;g++){
            if(g==0){
                for(k=0;k<q;k++)
                    for(j=0;j<q;j++)
                        res3[k*q+j] = theta[g][k*q+j]*n[g]/Psi[i+g*p];
            }else{
                for(k=0;k<q;k++)
                    for(j=0;j<q;j++)
                        res3[k*q+j] += theta[g][k*q+j]*n[g]/Psi[i+g*p];
            }
        }
        /* Invert theta sum */
        GaussJordan(q, res3, res1, det);

        /* Now solve for row i of lambda */
        for(j=0;j<q;j++)
        {result[j] = res2[i*q+j];}
        vec_mx_mult(q, q, result, res1, lambda0);
        for(j=0;j<q;j++)
        {lambda[i*q+j] = lambda0[j];}
    }

    free(tran); free(res1); free(res2); free(res3); free(result); free(lambda0);
}

/* For CUC */
void update_lambda_cuc(double *lambda, double **beta, double **s, double **theta, double *n, double *Psi, int p,
                       int q, int G){

    int i,j,g;
    double det[1];
    double *tran = malloc(sizeof(double)*p*q);
    double *res1 = malloc(sizeof(double)*p*q);
    double *res2 = malloc(sizeof(double)*p*q);
    double *res3 = malloc(sizeof(double)*q*q);


    /* sum_g[sampcov_g*beta_g'*n_g/psi_g] */
    for(g=0;g<G;g++){
        mx_trans(q,p,beta[g],tran);
        mx_mult(p,p,q,s[g],tran,res1);
        if(g==0){
            for(i=0;i<p;i++)
                for(j=0;j<q;j++)
                    res2[i*q+j] = res1[i*q+j]*n[g]/Psi[g];
        }else{
            for(i=0;i<p;i++)
                for(j=0;j<q;j++)
                    res2[i*q+j] += res1[i*q+j]*n[g]/Psi[g];
        }
    }

    /* sum_g[theta_g'*n_g/psi_g]^{-1} */
    for(g=0;g<G;g++){
        if(g==0){
            for(i=0;i<q;i++)
                for(j=0;j<q;j++)
                    res3[i*q+j] = theta[g][i*q+j]*n[g]/Psi[g];
        }else{
            for(i=0;i<q;i++)
                for(j=0;j<q;j++)
                    res3[i*q+j] += theta[g][i*q+j]*n[g]/Psi[g];
        }
    }

    GaussJordan(q, res3, res1,det); /* inverting the theta sum */
    mx_mult(p,q,q,res2,res1,lambda); /* this gives lambda */

    free(tran); free(res1); free(res2); free(res3);
}

// ---------- PSI UPDATES ---------- //

double update_psi(double *lambda, double *beta, double *sampcovtilde, int p, int q) {
    // performs the psi update for CCC & UUC models

    int i;
    double psi;
    double *result_1 = malloc(sizeof(double)*p*p);
    double *result_2 = malloc(sizeof(double)*p);

    // lambda * beta * sampcovtilde
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovtilde, result_2);

    psi = 0.0;
    for(i=0;i<p;i++)
        psi += sampcovtilde[i*p+i]-result_2[i];
    psi /= p;

    free(result_1); free(result_2);

    return psi;
}

void update_psi2(double *psi, double *lambda, double *beta, double *sampcovtilde, int p, int q){
    // updates psi for CCU & UUU models

    int i;
    double *result_1 = malloc(sizeof(double)*p*p);
    double *result_2 = malloc(sizeof(double)*p);

    // lambda * beta * sampcovtilde
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovtilde, result_2);

    for(i=0;i<p;i++){
        psi[i] = sampcovtilde[i*p+i]-result_2[i];
    }
    free(result_1); free(result_2);
}

double update_psi_cuc(double *lambda, double *beta, double *sampcovg, double *theta, int p, int q){

    int i;
    double psi;
    double *temp = malloc(sizeof(double)*q*p);
    double *result_1 = malloc(sizeof(double)*p*p);
    double *result_2 = malloc(sizeof(double)*p);
    double *result_3 = malloc(sizeof(double)*p);


    // lambda * beta * sampcov
    mx_mult(p,q,p,lambda, beta, result_1);
    mx_mult_diag1(p,p,result_1, sampcovg, result_2);

    // lambda * theta * lambda'
    mx_trans(p,q,lambda,temp);
    mx_mult(p,q,q,lambda,theta,result_1);
    mx_mult_diag1(p,q,result_1,temp,result_3);

    psi = 0.0;
    for(i=0;i<p;i++)
        psi += sampcovg[i*p+i]-2*result_2[i]+result_3[i];
    psi /= p;

    free(temp); free(result_1); free(result_2);
    free(result_3);

    return psi;
}

void update_psi_cuu(double *psi, double *lambda, double **beta, double **sampcovg, double **theta, int p, int q, int G){

    int i, g, j;
    double *temp = malloc(sizeof(double)*q*p);
    double *result_1 = malloc(sizeof(double)*p*p);
    double *result_2 = malloc(sizeof(double)*G*p);
    double *result_3 = malloc(sizeof(double)*G*p);
    double *result = malloc(sizeof(double)*p);


    // lambda * beta * sampcov
    for (g=0;g<G;g++){
        mx_mult(p,q,p,lambda, beta[g], result_1);
        mx_mult_diag1(p,p,result_1, sampcovg[g], result);
        for(j=0; j<p; j++) {result_2[g*p+j] = result[j];}
    }

    // lambda * theta * lambda'
    for(g=0;g<G;g++){
        mx_trans(p,q,lambda,temp);
        mx_mult(p,q,q,lambda,theta[g],result_1);
        mx_mult_diag1(p,q,result_1,temp,result);
        for(j=0; j<p; j++) {result_3[g*p+j] = result[j];}
    }

    for(g=0;g<G;g++)
        for(i=0;i<p;i++)
            psi[g*p+i] = sampcovg[g][i*p+i]-2*result_2[g*p+i]+result_3[g*p+i];

    free(temp); free(result_1); free(result_2);
    free(result_3); free(result);
}

void update_psi_ucu(double *psi, double **lambda, double **beta, double **sampcov, int p, int q, double *pi, int G){

    int i,g,j;
    double *result_1 = malloc(sizeof(double)*p*p);
    double *result_2 = malloc(sizeof(double)*G*p);
    double *result = malloc(sizeof(double)*p);


    /* lambda*beta*sampcovtilde */
    for (g=0; g<G; g++){
        mx_mult(p,q,p,lambda[g], beta[g], result_1);
        mx_mult_diag1(p,p,result_1, sampcov[g], result);
        for(j=0; j<p; j++) {result_2[g*p+j] = result[j];}
    }

    for(i=0;i<p;i++){
        psi[i]=0.0;
        for(g=0;g<G;g++) psi[i] += pi[g]*(sampcov[g][i*p+i]-result_2[i+g*p]);
    }

    free(result_1); free(result_2); free(result);
}

double update_psi_ucc(double **lambda, double **beta, double **sampcov, int p, int q, double *pi, int G){

    int i,g,j;
    double psi;
    double *result_1 = malloc(sizeof(double)*p*p);
    double *result_2 = malloc(sizeof(double)*G*p);
    double *result = malloc(sizeof(double)*p);

    // lambda * beta * sampcovtilde
    for (g=0; g<G; g++){
        mx_mult(p,q,p,lambda[g], beta[g], result_1);
        mx_mult_diag1(p,p,result_1, sampcov[g], result);
        for(j=0; j<p; j++) {result_2[g*p+j] = result[j];}
    }

    psi=0.0;
    for(g=0;g<G;g++)
        for(i=0;i<p;i++)
            psi += pi[g]*(sampcov[g][i*p+i]-result_2[g*p+i]);
    psi /= p;

    free(result_1); free(result_2); free(result);

    return psi;
}

// ---------- DETSIG UPDATES ---------- //

double update_detsig_iso(double *lambda, double psi, double log_detpsi, int p, int q){

    int i, j;
    double det[1];

    double *tmp = malloc(sizeof(double)*p*p);
    double *tmp2 = malloc(sizeof(double)*p*p);

    update_beta_iso(tmp2, psi, lambda, p, q);
    mx_mult(q,p,q,tmp2, lambda, tmp);
    for(i=0;i<q;i++){
        for(j=0;j<q;j++){
            tmp2[i*q+j] = tmp[i*q+j]*(-1.0);
            if(i==j) tmp2[i*q+i] += 1.0;
        }
    }
    GaussJordan(q,tmp2,tmp,det);

    free(tmp); free(tmp2);

    return (log_detpsi - log(det[0]));
}

double update_detsig_niso(double *lambda, double *psi, double log_detpsi, int p, int q){

    int i, j;
    double det[1];

    double *tmp = malloc(sizeof(double)*p*p);
    double *tmp2 = malloc(sizeof(double)*p*p);

    update_beta_niso(tmp2, psi, lambda, p, q);
    mx_mult(q,p,q,tmp2, lambda, tmp);
    for(i=0;i<q;i++){
        for(j=0;j<q;j++){
            tmp2[i*q+j] = tmp[i*q+j]*(-1.0);
            if(i==j) tmp2[i*q+i] += 1.0;
        }
    }
    GaussJordan(q,tmp2,tmp,det);

    free(tmp); free(tmp2);

    return (log_detpsi - log(det[0]));
}
