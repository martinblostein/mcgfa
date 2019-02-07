#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include "functions.h"

// AECM ALGORITHMS: One for each set of covariance constraints

double aecm_CCC(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N,
                  double *mu, double *lambda, double *psi_ptr, double *eta, double *alpha, int fixed_alpha, int fixed_eta,
                  double eta_max, double alpha_min, double tol, int max_it, int *iterations){
    int i, g;
    double bic;
    int it = 0, stop = 0, params;

    double psi, log_detpsi, log_c=0.0, log_detsig;
    double *pi = malloc(sizeof(double)*G);
    double *n = malloc(sizeof(double)*G);
    double *at = malloc(sizeof(double)*150000);
    double *l = malloc(sizeof(double)*150000);
    double *sampcovtilde = malloc(sizeof(double)*p*p);
    double *max_log_dens = malloc(sizeof(double)*N);
    double *log_dens = malloc(sizeof(double)*N*G);
    double *beta = malloc(sizeof(double)*q*p);
    double *theta = malloc(sizeof(double)*q*q);
    double *zbad = malloc(sizeof(double)*N*G);
    double *correction = malloc(sizeof(double)*N*G);
    double **mahal   = malloc(sizeof(double*)*G);
    for(g=0; g < G; g++) mahal[g] = malloc(sizeof(double)*N);

    psi = *psi_ptr;

    // initialize determinants & density constant
    log_detpsi = p*log(psi);
    log_detsig = update_detsig_iso(lambda, psi, log_detpsi, p, q);
    log_c = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig;

	// MAIN AECM LOOP
    while(stop==0){

        update_n(n, z, G, N);
        update_pi(pi, n, G, N);
        update_correction(correction, v, eta, G, N); // correction is (vig+(1-vig)/eta)
        update_mu(mu, n, x, z, correction, G, N, p);
        update_alpha(alpha, alpha_min, z, n, v, G, N, fixed_alpha);

        update_mahal_CCC(mahal, x, lambda, psi, mu, N, G, p, q);
        for(i=0;i<N;i++) for(g=0;g<G;g++) zbad[i*G+g] = z[i*G+g]*(1-v[i*G+g]);
        update_eta(eta, eta_max, zbad, mahal, N, G, p, fixed_eta);

        update_zv(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        update_correction(correction, v, eta, G, N);
        update_stilde(sampcovtilde,x, z, mu, correction, G, N, p); // sample cov update

        update_beta_iso(beta, psi, lambda, p, q);
        update_theta(theta, beta, lambda, sampcovtilde, p, q);
        update_lambda(lambda, beta, sampcovtilde, theta, p, q);
        psi = update_psi(lambda, beta, sampcovtilde, p, q);

        log_detpsi = p*log(psi);
        log_detsig = update_detsig_iso(lambda, psi, log_detpsi, p, q);
        log_c = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig;

        update_mahal_CCC(mahal, x, lambda, psi, mu, N, G, p, q);
        update_zv(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        // stopping criteria
        stop = converge_test(l, at, max_log_dens, log_dens, N, it, G, tol);
        it++;
        if (it >= max_it) {stop=1;}
    }

    *iterations=it;
    params = G-1 + G*p + p*q - q*(q-1)/2 + 1;
    params = params + (fixed_alpha ? 1 : G) + (fixed_eta ? 1 : G);

    bic = 2.0*l[it-1] - params*log(N);

    *psi_ptr = psi;

    free(n); free(beta); free(theta); free(sampcovtilde);
    free(l); free(at); free(pi); free(zbad); free(correction);
    free(max_log_dens); free(log_dens);
    
    for(g=0; g < G; g++) {
        free(mahal[g]);
    }
    
    free(mahal);

    return bic;
}

double aecm_CCU(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N,
                  double *mu, double *lambda, double *psi, double *eta, double *alpha, int fixed_eta, int fixed_alpha,
                  double eta_max, double alpha_min, double tol, int max_it, int *iterations) {
    int i, g;
    double bic;
    int it = 0, stop = 0, params;

    double log_detpsi, log_detsig, log_c=0.0;
    double *det = malloc(sizeof(double)*G);
    double *pi = malloc(sizeof(double)*G);
    double *n = malloc(sizeof(double)*G);
    double *at = malloc(sizeof(double)*150000);
    double *l = malloc(sizeof(double)*150000);
    double *sampcovtilde = malloc(sizeof(double)*p*p);
    double *beta = malloc(sizeof(double)*q*p);
    double *theta = malloc(sizeof(double)*q*q);
    double *w = malloc(sizeof(double)*G*N);
    double *max_log_dens = malloc(sizeof(double)*N);
    double *log_dens = malloc(sizeof(double)*N*G);
    double *zbad = malloc(sizeof(double)*N*G);
    double *correction = malloc(sizeof(double)*N*G);
    double **mahal   = malloc(sizeof(double*)*G);
    for(g=0; g < G; g++) mahal[g] = malloc(sizeof(double)*N);



    // initialize determinants & density constant
    log_detpsi = 0.0; for(i=0;i<p;i++) log_detpsi += log(psi[i]);
    log_detsig = update_detsig_niso(lambda, psi, log_detpsi, p, q);
    log_c = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig;

	// MAIN AECM LOOP
    while(stop==0){

        update_n(n, z, G, N);
        update_pi(pi, n, G, N);
        update_correction(correction, v, eta, G, N); // correction is (vig+(1-vig)/eta)
        update_mu(mu, n, x, z, correction, G, N, p);
        update_alpha(alpha, alpha_min, z, n, v, G, N, fixed_alpha);

        for(i=0;i<N;i++) for(g=0;g<G;g++) zbad[i*G+g] = z[i*G+g]*(1-v[i*G+g]);
        update_eta(eta, eta_max, zbad, mahal, N, G, p, fixed_eta);

        update_mahal_CCU(mahal, x, lambda, psi, mu, N, G, p, q);
        update_zv(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
		if (cls_ind) known_z(cls,z,N,G);

		update_correction(correction, v, eta, G, N);
        update_stilde(sampcovtilde, x, z, mu, correction, G, N, p);

        update_beta_niso(beta, psi, lambda, p, q);
        update_theta(theta, beta, lambda, sampcovtilde, p, q);
        update_lambda(lambda, beta, sampcovtilde, theta, p, q);
        update_psi2(psi, lambda, beta, sampcovtilde, p, q);

        log_detpsi = 0.0; for(i=0;i<p;i++) log_detpsi += log(psi[i]);
        log_detsig = update_detsig_niso(lambda, psi, log_detpsi, p, q);
        log_c = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig;

        update_mahal_CCU(mahal, x, lambda, psi, mu, N, G, p, q);
        update_zv(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        // stopping criteria
        stop = converge_test(l, at, max_log_dens, log_dens, N, it, G, tol);
        it++;
        if (it >= max_it) {stop=1;}
    }

    *iterations=it;

    params = G-1 + G*p + p*q - q*(q-1)/2 + p;
    params = params + (fixed_alpha ? 1 : G) + (fixed_eta ? 1 : G);

    bic = 2.0*l[it-1] - params*log(N);

    // Deallocate memory
    free(w);  free(n); free(det); free(beta);
    free(theta); free(sampcovtilde); free(l); free(at); free(pi);
    free(zbad); free(correction); free(max_log_dens);
    free(log_dens);

    for(g=0; g < G; g++) {
        free(mahal[g]);
    }
    
    free(mahal);
    
    return bic;
}

double aecm_CUC(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N,
                  double *mu, double *lambda, double *psi, double *eta, double *alpha, int fixed_eta, int fixed_alpha,
                  double eta_max, double alpha_min, double tol, int max_it, int *iterations){
    int i, g;
    double bic;
    int it = 0, stop = 0, params;

    double *log_c = malloc(sizeof(double)*G);
    double *log_detpsi = malloc(sizeof(double)*G);
    double *log_detsig = malloc(sizeof(double)*G);
    double *pi = malloc(sizeof(double)*G);
    double *n = malloc(sizeof(double)*G);
    double *at = malloc(sizeof(double)*150000);
    double *l = malloc(sizeof(double)*150000);
    double **sampcov    = malloc(sizeof(double*)*G);
    double **beta    = malloc(sizeof(double*)*G);
    double **theta    = malloc(sizeof(double*)*G);
    double **mahal   = malloc(sizeof(double*)*G);
    for(g=0; g < G; g++) {
        sampcov[g] = malloc(sizeof(double)*p*p);
        beta[g]    = malloc(sizeof(double)*q*p);
        theta[g]   = malloc(sizeof(double)*q*q);
        mahal[g]   = malloc(sizeof(double)*N);
    }
    double *max_log_dens = malloc(sizeof(double)*N);
    double *log_dens = malloc(sizeof(double)*N*G);
    double *zbad = malloc(sizeof(double)*N*G);
    double *correction = malloc(sizeof(double)*N*G);

    // initialize determinants & density constant
    for(g=0;g<G;g++){
        log_detpsi[g] = p*log(psi[g]);
        log_detsig[g] = update_detsig_iso(lambda, psi[g], log_detpsi[g], p, q);
        log_c[g] = (p/2.0)*log(2*M_PI) + 0.5*log_detsig[g];
    }

    // MAIN AECM LOOP
    while(stop==0){

        update_n(n, z, G, N);
        update_pi(pi, n, G, N);
        update_correction(correction, v, eta, G, N); // correction is (vig+(1-vig)/eta)
        update_mu(mu, n, x, z, correction, G, N, p);
        update_alpha(alpha, alpha_min, z, n, v, G, N, fixed_alpha);

        update_mahal_CUC(mahal, x, lambda, psi, mu, N, G, p, q);
        for(i=0;i<N;i++) for(g=0;g<G;g++) zbad[i*G+g] = z[i*G+g]*(1-v[i*G+g]);
        update_eta(eta, eta_max, zbad, mahal, N, G, p, fixed_eta);

        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        update_correction(correction, v, eta, G, N);
        update_n(n, z, G, N);
        update_sg(sampcov, x, z, mu, correction, n, p, G, N);

        for (g=0;g<G;g++) update_beta_iso(beta[g], psi[g], lambda, p, q);
        for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        update_lambda_cuc(lambda, beta, sampcov, theta, n, psi, p, q, G);
        for(g=0;g<G;g++) psi[g] = update_psi_cuc(lambda, beta[g], sampcov[g], theta[g], p, q);

        for(g=0;g<G;g++){
            log_detpsi[g] = p*log(psi[g]);
            log_detsig[g] = update_detsig_iso(lambda, psi[g], log_detpsi[g], p, q);
            log_c[g] = (p/2.0)*log(2*M_PI) + 0.5*log_detsig[g];
        }

        update_mahal_CUC(mahal, x, lambda, psi, mu, N, G, p, q);
        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        // stopping criteria
        stop = converge_test(l, at, max_log_dens, log_dens, N, it, G, tol);
        it++;
        if (it >= max_it) {stop=1;}
    }

    *iterations = it;

    params = G-1 + G*p+ p*q - q*(q-1)/2 + G;
    params = params + (fixed_alpha ? 1 : G) + (fixed_eta ? 1 : G);

    bic = 2.0*l[it-1] - params*log(N);

    /* Deallocate memory */

    free(log_dens); free(n); free(log_c); free(max_log_dens);
    free(l); free(at); free(pi); free(log_detpsi);
    free(zbad); free(correction);
    
    for(g=0; g < G; g++) {
        free(beta[g]);
        free(theta[g]);
        free(sampcov[g]);
        free(mahal[g]);
    }

    free(beta); free(theta); free(sampcov); free(mahal);

    return bic;
}

double aecm_CUU(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N,
                  double *mu, double *lambda, double *Psi, double *eta, double *alpha, int fixed_eta, int fixed_alpha,
                  double eta_max, double alpha_min, double tol, int max_it, int *iterations){

    int i, j, g;
    double bic;
    int it = 0, stop = 0, params;

    double *max_log_dens = malloc(sizeof(double)*N);
    double *log_dens = malloc(sizeof(double)*N*G);
    double *pi = malloc(sizeof(double)*G);
    double *n = malloc(sizeof(double)*G);
    double *at = malloc(sizeof(double)*150000);
    double *l = malloc(sizeof(double)*150000);
    double **sampcov    = malloc(sizeof(double*)*G);
    double **beta    = malloc(sizeof(double*)*G);
    double **theta    = malloc(sizeof(double*)*G);
    double **mahal   = malloc(sizeof(double*)*G);

    for(g=0; g < G; g++) {
        sampcov[g] = malloc(sizeof(double)*p*p);
        beta[g]    = malloc(sizeof(double)*q*p);
        theta[g]   = malloc(sizeof(double)*q*q);
        mahal[g]   = malloc(sizeof(double)*N);
    }
    double *w = malloc(sizeof(double)*G*N);
    double *log_detpsi = malloc(sizeof(double)*G);
    double *log_detsig = malloc(sizeof(double)*G);
    double *log_c = malloc(sizeof(double)*G);
    double *Psi0 = malloc(sizeof(double)*p);
    double *zbad = malloc(sizeof(double)*N*G);
    double *correction = malloc(sizeof(double)*N*G);

    // initialize determinants & density constant
    for (g=0;g<G;g++){
        log_detpsi[g]=0.0;
        for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g*p+j]);
    }
    for(g=0;g<G;g++){
        for(j=0; j<p; j++){
            Psi0[j] = Psi[g*p+j];
        }
        log_detsig[g] = update_detsig_niso(lambda, Psi0, log_detpsi[g], p, q);
    }
    for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

    // MAIN AECM LOOP
    while(stop==0){

        update_n(n, z, G, N);
        update_pi(pi, n, G, N);
        update_correction(correction, v, eta, G, N); // correction is (vig+(1-vig)/eta)
        update_mu(mu, n, x, z, correction, G, N, p);
        update_alpha(alpha, alpha_min, z, n, v, G, N, fixed_alpha);

        update_mahal_CUU(mahal, x, lambda, Psi, mu, N, G, p, q);
        for(i=0;i<N;i++) for(g=0;g<G;g++) zbad[i*G+g] = z[i*G+g]*(1-v[i*G+g]);
        update_eta(eta, eta_max, zbad, mahal, N, G, p, fixed_eta);

        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        update_n(n, z, G, N);
        update_correction(correction, v, eta, G, N);
        update_sg(sampcov, x, z, mu, correction, n, p, G, N);

        for(g=0;g<G;g++){
            for(j=0; j<p; j++){
                Psi0[j] = Psi[g*p+j];
            }

            update_beta_niso(beta[g], Psi0, lambda, p, q);
        }
        for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda, sampcov[g], p, q);
        update_lambda_cuu(lambda, beta, sampcov, theta, n, Psi, p, q, G);
        update_psi_cuu(Psi, lambda, beta, sampcov, theta, p, q, G);

        for (g=0;g<G;g++){
            log_detpsi[g]=0.0;
            for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g*p+j]);
        }
        for(g=0;g<G;g++){
            for(j=0; j<p; j++){
                Psi0[j] = Psi[g*p+j];
            }
            log_detsig[g] = update_detsig_niso(lambda, Psi0, log_detpsi[g], p, q);
        }
        for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_mahal_CUU(mahal, x, lambda, Psi, mu, N, G, p, q);
        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        // stopping criteria
        stop = converge_test(l, at, max_log_dens, log_dens, N, it, G, tol);
        it++;
        if (it >= max_it) {stop=1;}
    }

    *iterations = it;

    params = G-1 + G*p + p*q - q*(q-1)/2 + G*p;
    params = params + (fixed_alpha ? 1 : G) + (fixed_eta ? 1 : G);

    bic = 2.0*l[it-1] - params*log(N);

    /* Deallocate memory */
    free(w); free(n); free(l); free(at); free(pi);
    free(log_detsig); free(log_c); free(log_detpsi); free(Psi0);
    free(max_log_dens); free(log_dens); 
    free(zbad); free(correction);

    for(g=0; g < G; g++) {
        free(beta[g]);
        free(theta[g]);
        free(sampcov[g]);
        free(mahal[g]);
    }
    free(beta); free(theta); free(sampcov); free(mahal);
    
    return bic;
    }

double aecm_UCC(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N,
                  double *mu, double *lam_vec, double *psi_ptr, double *eta, double *alpha, int fixed_eta, int fixed_alpha,
                  double eta_max, double alpha_min, double tol, int max_it, int *iterations){
    int i, g;
    double bic;
    int it = 0, stop = 0, params;

    double log_detpsi, psi;
    double *max_log_dens = malloc(sizeof(double)*N);
    double *log_dens = malloc(sizeof(double)*N*G);
    double *log_c = malloc(sizeof(double)*G);
    double *log_detsig = malloc(sizeof(double)*G);
    double *pi = malloc(sizeof(double)*G);
    double *n = malloc(sizeof(double)*G);
    double *at = malloc(sizeof(double)*150000);
    double *l = malloc(sizeof(double)*150000);
    double **sampcov    = malloc(sizeof(double*)*G);
    double **lambda    = malloc(sizeof(double*)*G);
    double **beta    = malloc(sizeof(double*)*G);
    double **theta    = malloc(sizeof(double*)*G);
    double **mahal   = malloc(sizeof(double*)*G);
    for(g=0; g < G; g++) {
        sampcov[g] = malloc(sizeof(double)*p*p);
        lambda[g]  = malloc(sizeof(double)*p*q);
        beta[g]    = malloc(sizeof(double)*q*p);
        theta[g]   = malloc(sizeof(double)*q*q);
        mahal[g]   = malloc(sizeof(double)*N);
    }

    double *w = malloc(sizeof(double)*G*N);
    double *zbad = malloc(sizeof(double)*N*G);
    double *correction = malloc(sizeof(double)*N*G);

    psi = *psi_ptr;
    get_lambda(lam_vec,lambda,G,p,q);

    // initialize determinants & density constant
    log_detpsi = 0.0;
    for(i=0;i<p;i++) log_detpsi += log(psi);
    for(g=0;g<G;g++) log_detsig[g] = update_detsig_iso(lambda[g], psi, log_detpsi, p, q);
    for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

    while(stop==0){

        update_n(n, z, G, N);
        update_pi(pi, n, G, N);
        update_correction(correction, v, eta, G, N); // correction is (vig+(1-vig)/eta)
        update_mu(mu, n, x, z, correction, G, N, p);
        update_alpha(alpha, alpha_min, z, n, v, G, N, fixed_alpha);

        update_mahal_UCC(mahal, x, lambda, psi, mu, N, G, p, q);
        for(i=0;i<N;i++) for(g=0;g<G;g++) zbad[i*G+g] = z[i*G+g]*(1-v[i*G+g]);
        update_eta(eta, eta_max, zbad, mahal, N, G, p, fixed_eta);

        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        update_n(n, z, G, N);
        update_correction(correction, v, eta, G, N);
        update_sg(sampcov, x, z, mu, correction, n, p, G, N);

        for(g=0;g<G;g++) update_beta_iso(beta[g], psi, lambda[g], p, q);
        for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);

        psi = update_psi_ucc(lambda, beta, sampcov, p, q, pi, G);

        log_detpsi = 0.0;
        for(i=0;i<p;i++) log_detpsi += log(psi);
        for(g=0;g<G;g++) log_detsig[g] = update_detsig_iso(lambda[g], psi, log_detpsi, p, q);
        for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_mahal_UCC(mahal, x, lambda, psi, mu, N, G, p, q);
        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        // stopping criteria
        stop = converge_test(l, at, max_log_dens, log_dens, N, it, G, tol);
        it++;
        if (it >= max_it) {stop=1;}
    }

    *iterations = it;
    params = G-1 + G*p + G*(p*q - q*(q-1)/2) + 1;
    params = params + (fixed_alpha ? 1 : G) + (fixed_eta ? 1 : G);
    
    bic = 2.0*l[it-1] - params*log(N);

    store_lambda(lam_vec, lambda, G, p, q);

    *psi_ptr = psi;
    free(w); free(n); free(l); free(at); free(pi); free(log_detsig); free(log_c);
    free(zbad); free(correction); free(max_log_dens); free(log_dens);
    
    for(g=0; g < G; g++) {
        free(beta[g]);
        free(lambda[g]);
        free(theta[g]);
        free(sampcov[g]);
        free(mahal[g]);
    }

    free(beta); free(lambda); free(theta); free(sampcov); free(mahal);

    return bic;
}

double aecm_UCU(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N,
                  double *mu, double *lam_vec, double *psi, double *eta, double *alpha, int fixed_eta, int fixed_alpha,
                  double eta_max, double alpha_min, double tol, int max_it, int *iterations){
    int i, j=0, g;
    double bic;
    int it = 0, stop = 0, params;

    double log_detpsi;
    double *max_log_dens = malloc(sizeof(double)*N);
    double *log_dens = malloc(sizeof(double)*N*G);
    double *log_detsig = malloc(sizeof(double)*G);
    double *log_c = malloc(sizeof(double)*G);
    double *pi = malloc(sizeof(double)*G);
    double *n = malloc(sizeof(double)*G);
    double *at = malloc(sizeof(double)*150000);
    double *l = malloc(sizeof(double)*150000);
    double **sampcov    = malloc(sizeof(double*)*G);
    double **lambda    = malloc(sizeof(double*)*G);
    double **beta    = malloc(sizeof(double*)*G);
    double **theta    = malloc(sizeof(double*)*G);
    double **mahal   = malloc(sizeof(double*)*G);
    for(g=0; g < G; g++) {
        sampcov[g] = malloc(sizeof(double)*p*p);
        lambda[g]  = malloc(sizeof(double)*p*q);
        beta[g]    = malloc(sizeof(double)*q*p);
        theta[g]   = malloc(sizeof(double)*q*q);
        mahal[g]   = malloc(sizeof(double)*N);
    }
    double *zbad = malloc(sizeof(double)*N*G);
    double *correction = malloc(sizeof(double)*N*G);

    get_lambda(lam_vec,lambda,G,p,q);

    // initialize determinants & density constant
    log_detpsi = 0.0;
    for(j=0;j<p;j++) log_detpsi += log(psi[j]);
    for(g=0;g<G;g++) log_detsig[g] = update_detsig_niso(lambda[g], psi, log_detpsi, p, q);
    for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

    while(stop==0){

        update_n(n, z, G, N);
        update_pi(pi, n, G, N);
        update_correction(correction, v, eta, G, N); // correction is (vig+(1-vig)/eta)
        update_mu(mu, n, x, z, correction, G, N, p);
        update_alpha(alpha, alpha_min, z, n, v, G, N, fixed_alpha);

        update_mahal_UCU(mahal, x, lambda, psi, mu, N, G, p, q);
        for(i=0;i<N;i++) for(g=0;g<G;g++) zbad[i*G+g] = z[i*G+g]*(1-v[i*G+g]);
        update_eta(eta, eta_max, zbad, mahal, N, G, p, fixed_eta);

        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        update_n(n, z, G, N);
        update_correction(correction, v, eta, G, N);
        update_sg(sampcov, x, z, mu, correction, n, p, G, N);

        for(g=0;g<G;g++) update_beta_niso(beta[g], psi, lambda[g], p, q);
        for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);

        update_psi_ucu(psi, lambda, beta, sampcov, p, q, pi, G);

        log_detpsi = 0.0;
        for(j=0;j<p;j++) log_detpsi += log(psi[j]);
        for(g=0;g<G;g++) log_detsig[g] = update_detsig_niso(lambda[g], psi, log_detpsi, p, q);
        for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_mahal_UCU(mahal, x, lambda, psi, mu, N, G, p, q);
        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        // stopping criteria
        stop = converge_test(l, at, max_log_dens, log_dens, N, it, G, tol);
        it++;
        if (it >= max_it) {stop=1;}
    }

    *iterations = it;
    params = G-1 + G*p + G*(p*q - q*(q-1)/2) + p;
    params = params + (fixed_alpha ? 1 : G) + (fixed_eta ? 1 : G);
    
    bic = 2.0*l[it-1] - params*log(N);

    store_lambda(lam_vec, lambda, G, p, q);

    free(log_dens); free(n); free(max_log_dens); free(l); free(at);
    free(pi); free(log_detsig); free(log_c); free(zbad); free(correction);
    
    for(g=0; g < G; g++) {
        free(beta[g]);
        free(lambda[g]);
        free(theta[g]);
        free(sampcov[g]);
        free(mahal[g]);
    }
    free(beta); free(lambda); free(theta); free(sampcov); free(mahal);

    return bic;
}

double aecm_UUC(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N,
                  double *mu, double *lam_vec, double *psi, double *eta, double *alpha, int fixed_eta, int fixed_alpha,
                  double eta_max, double alpha_min, double tol, int max_it, int *iterations) {
    int i, g;
    double bic;
    int it = 0, stop = 0, params;

    double *max_log_dens = malloc(sizeof(double)*N);
    double *log_dens = malloc(sizeof(double)*N*G);
    double *log_detpsi = malloc(sizeof(double)*G);
    double *log_detsig = malloc(sizeof(double)*G);
    double *log_c = malloc(sizeof(double)*G);
    double *pi = malloc(sizeof(double)*G);
    double *n = malloc(sizeof(double)*G);
    double *at = malloc(sizeof(double)*150000);
    double *l = malloc(sizeof(double)*150000);
    double **sampcov    = malloc(sizeof(double*)*G);
    double **lambda    = malloc(sizeof(double*)*G);
    double **beta    = malloc(sizeof(double*)*G);
    double **theta    = malloc(sizeof(double*)*G);
    double **mahal   = malloc(sizeof(double*)*G);
    for(g=0; g < G; g++) {
        sampcov[g] = malloc(sizeof(double)*p*p);
        lambda[g]  = malloc(sizeof(double)*p*q);
        beta[g]    = malloc(sizeof(double)*q*p);
        theta[g]   = malloc(sizeof(double)*q*q);
        mahal[g]   = malloc(sizeof(double)*N);
    }
    double *zbad = malloc(sizeof(double)*N*G);
    double *correction = malloc(sizeof(double)*N*G);

    get_lambda(lam_vec,lambda,G,p,q);

    // initialize determinants & density constant
    for (g=0;g<G;g++) log_detpsi[g] = p*log(psi[g]);
    for (g=0;g<G;g++) log_detsig[g] = update_detsig_iso(lambda[g], psi[g], log_detpsi[g], p, q);
    for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

    // MAIN AECM LOOP
    while(stop==0) {

        update_n(n, z, G, N);
        update_pi(pi, n, G, N);
        update_correction(correction, v, eta, G, N); // correction is (vig+(1-vig)/eta)
        update_mu(mu, n, x, z, correction, G, N, p);
        update_alpha(alpha, alpha_min, z, n, v, G, N, fixed_alpha);

        update_mahal_UUC(mahal, x, lambda, psi, mu, N, G, p, q);
        for(i=0;i<N;i++) for(g=0;g<G;g++) zbad[i*G+g] = z[i*G+g]*(1-v[i*G+g]);
        update_eta(eta, eta_max, zbad, mahal, N, G, p, fixed_eta);

        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        update_n(n, z, G, N);
        update_correction(correction, v, eta, G, N);
        update_sg(sampcov, x, z, mu, correction, n, p, G, N);

        for(g=0;g<G;g++) update_beta_iso(beta[g], psi[g], lambda[g], p, q);
        for(g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
        for (g=0;g<G;g++) psi[g] = update_psi(lambda[g], beta[g], sampcov[g], p, q);

        for(g=0;g<G;g++) log_detpsi[g] = p*log(psi[g]);
        for(g=0;g<G;g++) log_detsig[g] = update_detsig_iso(lambda[g], psi[g], log_detpsi[g], p, q);
        for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_mahal_UUC(mahal, x, lambda, psi, mu, N, G, p, q);
        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        // stopping criteria
        stop = converge_test(l, at, max_log_dens, log_dens, N, it, G, tol);
        it++;
        if (it >= max_it) {stop=1;}
    }

    *iterations = it;
    params = G-1 + G*p + G*(p*q - q*(q-1)/2) + G;
    params = params + (fixed_alpha ? 1 : G) + (fixed_eta ? 1 : G);
    
    bic = 2.0*l[it-1] - params*log(N);

    store_lambda(lam_vec, lambda, G, p, q);

    free(log_dens); free(n); free(l); free(at); free(pi);  free(log_detpsi);
    free(log_detsig); free(log_c); free(zbad); free(correction);
    free(max_log_dens);
    
    for(g=0; g < G; g++) {
        free(beta[g]);
        free(lambda[g]);
        free(theta[g]);
        free(sampcov[g]);
        free(mahal[g]);
    }

    // memory deallocation
    free(beta); free(lambda); free(theta); free(sampcov); free(mahal);
    return bic;
}

double aecm_UUU(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N,
                  double *mu, double *lam_vec, double *Psi, double *eta, double *alpha, int fixed_eta, int fixed_alpha,
                  double eta_max, double alpha_min, double tol, int max_it, int *iterations) {

    int i, j=0, g;
    double bic;
    int it = 0, stop = 0, params;

    double *max_log_dens = malloc(sizeof(double)*N);
    double *log_dens = malloc(sizeof(double)*N*G);
    double *pi = malloc(sizeof(double)*G);
    double *n = malloc(sizeof(double)*G);
    double *at = malloc(sizeof(double)*150000);
    double *l = malloc(sizeof(double)*150000);
    double **sampcov = malloc(sizeof(double*)*G);
    double **lambda  = malloc(sizeof(double*)*G);
    double **beta    = malloc(sizeof(double*)*G);
    double **theta   = malloc(sizeof(double*)*G);
    double **mahal   = malloc(sizeof(double*)*G);
    for(g=0; g < G; g++) {
        sampcov[g] = malloc(sizeof(double)*p*p);
        lambda[g]  = malloc(sizeof(double)*p*q);
        beta[g]    = malloc(sizeof(double)*q*p);
        theta[g]   = malloc(sizeof(double)*q*q);
        mahal[g]   = malloc(sizeof(double)*N);
    }
    double *log_detpsi = malloc(sizeof(double)*G);
    double *log_detsig = malloc(sizeof(double)*G);
    double *log_c = malloc(sizeof(double)*G);
    double *Psi0 = malloc(sizeof(double)*p);
    double *zbad = malloc(sizeof(double)*N*G);
    double *correction = malloc(sizeof(double)*N*G);


    get_lambda(lam_vec,lambda,G,p,q);

    // initialize determinants & density constant
    for(g=0;g<G;g++) {
        log_detpsi[g]=0.0;
        for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g*p+j]);
    }
    for(g=0;g<G;g++) {
        for(j=0; j<p; j++) {
            Psi0[j] = Psi[g*p+j];
        }
        log_detsig[g] = update_detsig_niso(lambda[g], Psi0, log_detpsi[g], p, q);
    }
    for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];


    // MAIN AECM LOOP

    while(stop==0) {

        update_n(n, z, G, N);
        update_pi(pi, n, G, N);
        update_correction(correction, v, eta, G, N); // correction is (vig+(1-vig)/eta)
        update_mu(mu, n, x, z, correction, G, N, p);
        update_alpha(alpha, alpha_min, z, n, v, G, N, fixed_alpha);

        update_mahal_UUU(mahal, x, lambda, Psi, mu, N, G, p, q);
        for(i=0;i<N;i++) {for(g=0;g<G;g++) {zbad[i*G+g] = z[i*G+g]*(1-v[i*G+g]);}}
        update_eta(eta, eta_max, zbad, mahal, N, G, p, fixed_eta);

        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        update_correction(correction, v, eta, G, N);
        update_n(n,z,G,N);
        update_sg(sampcov, x, z, mu, correction, n, p, G, N);

        for(g=0;g<G;g++) {
            for(j=0; j<p; j++){
                Psi0[j] = Psi[g*p+j];
            }
            update_beta_niso(beta[g], Psi0, lambda[g], p, q);
        }

        for (g=0;g<G;g++) update_theta(theta[g], beta[g], lambda[g], sampcov[g], p, q);
        for (g=0;g<G;g++) update_lambda(lambda[g], beta[g], sampcov[g], theta[g], p, q);
        for (g=0;g<G;g++){
            update_psi2(Psi0, lambda[g], beta[g], sampcov[g], p, q);
            for(j=0; j<p; j++){
                Psi[g*p+j] = Psi0[j];
            }
        }

        for(g=0;g<G;g++) {
            log_detpsi[g]=0.0;
            for(j=0;j<p;j++) log_detpsi[g] += log(Psi[g*p+j]);
        }
        for(g=0;g<G;g++) {
            for(j=0; j<p; j++) {
                Psi0[j] = Psi[g*p+j];
            }
            log_detsig[g] = update_detsig_niso(lambda[g], Psi0, log_detpsi[g], p, q);

        }
        for (g=0;g<G;g++) log_c[g] = (p/2.0)*log(2.0*M_PI) + 0.5*log_detsig[g];

        update_mahal_UUU(mahal, x, lambda, Psi, mu, N, G, p, q);
        update_zv2(log_dens, x, z, v, pi, max_log_dens, log_c, mahal, eta, alpha, N, G, p, q);
        if (cls_ind) known_z(cls,z,N,G);

        // stopping criteria
        stop = converge_test(l, at, max_log_dens, log_dens, N, it, G, tol);
        it++;
        if (it >= max_it) stop =1 ;
    }

    *iterations = it;
    params = G-1 + G*p + G*(p*q - q*(q-1)/2) + G*p;
    params = params + (fixed_alpha ? 1 : G) + (fixed_eta ? 1 : G);
    
    bic = 2.0*l[it-1] - params*log(N);

    // copy parameter data back out of aecm function
    store_lambda(lam_vec, lambda, G, p, q);

    // memory deallocation
    free(log_dens); free(n); free(log_detpsi); free(l); free(at); free(pi);
    free(log_detsig); free(log_c); free(max_log_dens); free(Psi0);
    free(zbad); free(correction);

    for(g=0; g < G; g++) {
        free(beta[g]);
        free(lambda[g]);
        free(theta[g]);
        free(sampcov[g]);
        free(mahal[g]);
    }
    free(beta); free(lambda); free(theta); free(sampcov); free(mahal);
    return bic;
}

// Aitken criterion for convergence

int converge_test(double *l, double *at, double *v_max, double *log_dens, int N, int it, int G, double TOL){
    int i,g, flag=0;
    double sum, l_inf;

    l[it]=0.0;
    for(i=0; i<N; i++){
        sum=0.0;
        for(g=0; g<G; g++){
            sum += exp(log_dens[i*G+g]-v_max[i]);
        }
        l[it] += log(sum)+v_max[i];
        if(isnan(l[it])||isinf(l[it])) return -1;
    }

    if(it > 0)
        if(l[it]<l[it-1])
            return -1;

        if(it > 2){
            at[it-1]=(l[it]-l[it-1])/(l[it-1]-l[it-2]);
            if(at[it-1]<1.0){
                l_inf = l[it-1]+(l[it]-l[it-1])/(1-at[it-1]);
                if(fabs(l_inf - l[it])<TOL) flag=1;
            }
        }
        return flag;
}

// Functions to copy between 1D and 2D arrays

void get_lambda(double *lam_vec, double **lambda, int G, int p, int q) {
    // lambda is passed in as a 1D array. For the algorithms with unconstrained
    // lambda (beginning in U), lambda is copied into a 2D array
    int g,i,k=0;
    for (g=0;g<G;g++){
        for (i=0;i<p*q;i++){
            lambda[g][i] = lam_vec[k];
            k++;
        }
    }
}

void store_lambda(double *lam_vec, double **lambda, int G, int p, int q) {
    // this copies the 2D array back into the original 1D array
    int i,g,k=0;
    for(g=0;g<G;g++){
        for (i=0;i<p*q;i++){
            lam_vec[k]=lambda[g][i];
            k++;
        }
    }
}

// for labelled observations

void known_z(int *class, double *z, int N, int G) {
    // restores given labels after z-update
    int i,g;
    for(i=0;i<N;i++){
        if(class[i]!=0){
            for(g=1;g<=G;g++){
                z[i*G+g-1]=0.0;
                if(g==class[i]) z[i*G+g-1]=1.0;
            }
        }
    }
}
