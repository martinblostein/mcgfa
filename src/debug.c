#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"
#define OUT_DIR "/home/martin/Desktop/mcnfa_output/"

void clear_file(char *name) {
    char rem_file[80];
    strcpy(rem_file, OUT_DIR);
    strcat(rem_file, name);
    strcat(rem_file,".txt");

    FILE *fp;
    fp = fopen(rem_file,"w");
    fclose(fp);
}

void clear_all_files() {
    clear_file("mu");
    clear_file("alpha");
    clear_file("eta");
    clear_file("v");
    clear_file("z");
    clear_file("lambda");
    clear_file("Sigma");
    clear_file("correction");
    clear_file("snapshot");
    clear_file("BEFORE");
    clear_file("AFTER");
    clear_file("log_density");
    clear_file("S");
    clear_file("zbad");
}

void write_2D(double *arr, int n, int m, char *name, char *message, int replace) {
    // Prints n-by-m array stored in row-major order to output file with given name.
    char out_file[80];
    strcpy(out_file, OUT_DIR);
    strcat(out_file, name);
    strcat(out_file,".txt");

    FILE *fp;
    fp = (replace) ? fopen(out_file,"w") : fopen(out_file,"a");

    if (message[0] != '0') fprintf(fp,"%s\n\n",message);

    int i, j;
    for (i=0;i<n;i++) {
        for (j=0;j<m;j++) {
            fprintf(fp,"%.8f\t",arr[i*m+j]);
        }
        fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
    fclose(fp);
}

void write_3D(double **arr, int n, int m, int G, char *name, char *message) {
    int g;

    for (g=0; g<G; g++) {
        write_2D(arr[g], n, m, name, (g==0) ? message : "0", 0);
    }
}

void write_init(double *z, double *x, double *v, double **lambda, double *alpha, double *eta,
                int p, int q, int G, int N) {
    // prints starting values
    char *filename="initial";
    write_2D(x, N, p, filename, "initial X:", 1);
    write_2D(z, N, G, filename, "initial z:", 0);
    write_2D(v, N, G, filename, "initial v:", 0);
    write_2D(alpha, 1, G, filename, "initial alpha:", 0);
    write_2D(eta, 1, G, filename, "initial eta:", 0);
    write_3D(lambda, p, q, G, filename, "initial lambda");
}

void write_init_ccc(double *z, double *x, double *v, double *lambda, double *psi, double *alpha, double *eta,
                int p, int q, int G, int N) {
    // prints starting values
    char *filename="initial";
    write_2D(x, N, p, filename, "initial X:", 1);
    write_2D(z, N, G, filename, "initial z:", 0);
    write_2D(v, N, G, filename, "initial v:", 0);
    write_2D(alpha, 1, G, filename, "initial alpha:", 0);
    write_2D(eta, 1, G, filename, "initial eta:", 0);
    write_2D(lambda, p, q, filename, "initial lambda", 0);
    write_2D(psi,1,1, filename, "psi:", 0);
}

void write_init_ccu(double *z, double *x, double *v, double *lambda, double *psi, double *alpha, double *eta,
                    int p, int q, int G, int N) {
    // prints starting values
    char *filename="initial";
    write_2D(x, N, p, filename, "initial X:", 1);
    write_2D(z, N, G, filename, "initial z:", 0);
    write_2D(v, N, G, filename, "initial v:", 0);
    write_2D(alpha, 1, G, filename, "initial alpha:", 0);
    write_2D(eta, 1, G, filename, "initial eta:", 0);
    write_2D(lambda, p, q, filename, "initial lambda", 0);
    write_2D(psi,1,p, filename, "psi:", 0);
}


void snapshot(char *filename, int it, double *z, double *v, double *correction, double *mu, double *pi, double **lambda, double *alpha, double *eta,
              double *psi, double **sampcov, double **beta, double **theta, double *log_dens, int p, int q, int G, int N) {
    char it_message[80];
    sprintf(it_message,"ITERATION %d",it);
    write_2D(z, 0, 0, filename, it_message,0);
    write_2D(z, N, G, filename, "z:", 0);
    write_2D(v, N, G, filename, "v:", 0);
    write_2D(correction, N, G, filename, "correction:", 0);
    write_2D(alpha, 1, G, filename, "alpha:", 0);
    write_2D(eta, 1, G, filename, "eta:", 0);
    write_2D(mu, G, p, filename, "mu:", 0);
    write_2D(pi, 1, G, filename, "pi:", 0);
    write_3D(lambda, p, q, G, filename, "lambda:");
    write_2D(psi, G, p, filename, "psi: ", 0);
    write_3D(sampcov, p, p, G, filename, "S:");
    write_3D(beta, q, p, G, filename, "beta:");
    write_3D(theta, q, q, G, filename, "theta:");
    write_2D(log_dens, N, G, filename, "log_density: ", 0);
}

void snapshot_ccc(char *filename, int it, double *z, double *v, double *correction, double *mu, double *pi, double *lambda, double *alpha, double *eta,
              double *psi, double *sampcov, double *beta, double *theta, double *log_dens, int p, int q, int G, int N) {
    char it_message[80];
    sprintf(it_message,"ITERATION %d",it);
    write_2D(z, 0, 0, filename, it_message,0);
    write_2D(z, N, G, filename, "z:", 0);
    write_2D(v, N, G, filename, "v:", 0);
    write_2D(correction, N, G, filename, "correction:", 0);
    write_2D(alpha, 1, G, filename, "alpha:", 0);
    write_2D(eta, 1, G, filename, "eta:", 0);
    write_2D(mu, G, p, filename, "mu:", 0);
    write_2D(pi, 1, G, filename, "pi:", 0);
    write_2D(lambda, p, q, filename, "lambda:", 0);
    write_2D(psi, 1, 1, filename, "psi: ", 0);
    write_2D(sampcov, p, p, filename, "S:", 0);
    write_2D(beta, q, p, filename, "beta:", 0);
    write_2D(theta, q, q, filename, "theta:", 0);
    write_2D(log_dens, N, G, filename, "log_density: ", 0);
}

void snapshot_ccu(char *filename, int it, double *z, double *v, double *correction, double *mu, double *pi, double *lambda, double *alpha, double *eta,
                  double *psi, double *sampcov, double *beta, double *theta, double *log_dens, int p, int q, int G, int N) {
    char it_message[80];
    sprintf(it_message,"ITERATION %d",it);
    write_2D(z, 0, 0, filename, it_message,0);
    write_2D(z, N, G, filename, "z:", 0);
    write_2D(v, N, G, filename, "v:", 0);
    write_2D(correction, N, G, filename, "correction:", 0);
    write_2D(alpha, 1, G, filename, "alpha:", 0);
    write_2D(eta, 1, G, filename, "eta:", 0);
    write_2D(mu, G, p, filename, "mu:", 0);
    write_2D(pi, 1, G, filename, "pi:", 0);
    write_2D(lambda, p, q, filename, "lambda:", 0);
    write_2D(psi, 1, p, filename, "psi: ", 0);
    write_2D(sampcov, p, p, filename, "S:", 0);
    write_2D(beta, q, p, filename, "beta:", 0);
    write_2D(theta, q, q, filename, "theta:", 0);
    write_2D(log_dens, N, G, filename, "log_density: ", 0);
}

void snapshot_uuc(char *filename, int it, double *z, double *v, double *correction, double *mu, double *pi, double **lambda, double *alpha, double *eta,
                  double *psi, double **sampcov, double **beta, double **theta, double *log_dens, int p, int q, int G, int N) {
    char it_message[80];
    sprintf(it_message,"ITERATION %d",it);
    write_2D(z, 0, 0, filename, it_message,0);
    write_2D(z, N, G, filename, "z:", 0);
    write_2D(v, N, G, filename, "v:", 0);
    write_2D(correction, N, G, filename, "correction:", 0);
    write_2D(alpha, 1, G, filename, "alpha:", 0);
    write_2D(eta, 1, G, filename, "eta:", 0);
    write_2D(mu, G, p, filename, "mu:", 0);
    write_2D(pi, 1, G, filename, "pi:", 0);
    write_3D(lambda, p, q, G, filename, "lambda:");
    write_2D(psi, 1, G, filename, "psi: ", 0);
    write_3D(sampcov, p, p, G, filename, "S:");
    write_3D(beta, q, p, G, filename, "beta:");
    write_3D(theta, q, q, G, filename, "theta:");
    write_2D(log_dens, N, G, filename, "log_density: ", 0);
}
