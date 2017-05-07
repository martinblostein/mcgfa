// AECM algorithms

double aecm_CCC(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N, double *mu,
                  double *lam_vec, double *psi_vec, double *eta, double *alpha, double eta_max, double alpha_min,
                  double tol, int max_it, int *iterations);
double aecm_CCU(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N, double *mu,
                  double *lam_vec, double *psi_vec, double *eta, double *alpha, double eta_max, double alpha_min,
                  double tol, int max_it, int *iterations);
double aecm_CUC(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N, double *mu,
                  double *lam_vec, double *psi_vec, double *eta, double *alpha, double eta_max, double alpha_min,
                  double tol, int max_it, int *iterations);
double aecm_CUU(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N, double *mu,
                  double *lam_vec, double *psi_vec, double *eta, double *alpha, double eta_max, double alpha_min,
                  double tol, int max_it, int *iterations);

double aecm_UCC(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N, double *mu,
                  double *lam_vec, double *psi_vec, double *eta, double *alpha, double eta_max, double alpha_min,
                  double tol, int max_it, int *iterations);
double aecm_UCU(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N, double *mu,
                  double *lam_vec, double *psi_vec, double *eta, double *alpha, double eta_max, double alpha_min,
                  double tol, int max_it, int *iterations);
double aecm_UUC(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N, double *mu,
                  double *lam_vec, double *psi_vec, double *eta, double *alpha, double eta_max, double alpha_min,
                  double tol, int max_it, int *iterations);
double aecm_UUU(double *z, double *x, double *v, int cls_ind, int* cls, int q, int p, int G, int N, double *mu,
                  double *lam_vec, double *psi_vec, double *eta, double *alpha, double eta_max, double alpha_min,
                  double tol, int max_it, int *iterations);

// z-updates

void update_zv(double *log_dens, double *x, double *z, double *v, double *pi, double *max_log_dens, double log_c, double **mahal,
               double *eta, double *alpha, int N, int G, int p, int q);
void update_zv2(double *log_dens, double *x, double *z, double *v, double *pi, double *max_log_dens, double *log_c, double **mahal,
               double *eta, double *alpha, int N, int G, int p, int q);

// mahal-updates (basically woodbury output used in both zv- and eta- updates)

double Brent_fmin(double ax, double bx, double (*f)(double, void *), void *info, double tol);

void update_eta(double *eta, double eta_max, double *zbad, double **mahalanobis, int N, int G, int p);

void update_mahal_CCC(double **mahal, double *x, double  *lambda, double  psi, double *mu, int N, int G, int p, int q);
void update_mahal_CCU(double **mahal, double *x, double  *lambda, double *psi, double *mu, int N, int G, int p, int q);
void update_mahal_CUC(double **mahal, double *x, double  *lambda, double *psi, double *mu, int N, int G, int p, int q);
void update_mahal_CUU(double **mahal, double *x, double  *lambda, double *psi, double *mu, int N, int G, int p, int q);

void update_mahal_UCC(double **mahal, double *x, double **lambda, double  psi, double *mu, int N, int G, int p, int q);
void update_mahal_UCU(double **mahal, double *x, double **lambda, double *psi, double *mu, int N, int G, int p, int q);
void update_mahal_UUC(double **mahal, double *x, double **lambda, double *psi, double *mu, int N, int G, int p, int q);
void update_mahal_UUU(double **mahal, double *x, double **lambda, double *psi, double *mu, int N, int G, int p, int q);

// other updates

void   update_n(double *n, double *z, int G, int N);
void   update_pi(double *pi, double *n, int G, int N);
void   update_mu(double *mu, double *n, double *x, double *z, double *correction, int G, int N, int p);
void   update_correction(double *correction, double *v, double *eta, int G, int N);

void   update_alpha(double *alpha, double alpha_min, double *z, double *n, double *v, int G, int N);
void   update_alpha_numerical(double *alpha, double alpha_min, double *z, double *v, int G, int N);

void   update_stilde(double *sampcovtilde, double *x, double *z, double *mu, double *correction, int G, int N, int p);
void   update_sg(double **S, double *x, double *z, double *mu, double *correction, double *n, int p, int G, int N);

void   update_beta_iso(double *beta, double psi, double *lambda, int p, int q);
void   update_beta_niso(double *beta, double *PSI, double *lambda, int p, int q);

void   update_theta(double *theta, double *beta, double *lambda, double *sampcovtilde, int p, int q);

void   update_lambda(double *lambda, double *beta, double *s, double *theta, int p, int q);
void   update_lambda_cuc(double *lambda, double **beta, double **s, double **theta, double *n, double *Psi, int p, int q, int G);
void   update_lambda_cuu(double *lambda, double **beta, double **s, double **theta, double *n, double *Psi, int p, int q, int G);

double update_psi(double *lambda, double *beta, double *sampcovtilde, int p, int q);
void   update_psi2(double *psi, double *lambda, double *beta, double *sampcovtilde, int p, int q);
double update_psi_cuc(double *lambda, double *beta, double *sampcovg, double *theta, int p, int q);
void   update_psi_ucu(double *psi, double **lambda, double **beta, double **sampcov, int p, int q, double *pi, int G);
double update_psi_ucc(double **lambda, double **beta, double **sampcov, int p, int q, double *pi, int G);
void   update_psi_cuu(double *psi, double *lambda, double **beta, double **sampcovg, double **theta, int p, int q, int G);

double update_detsig_iso(double *lambda, double psi, double log_detpsi, int p, int q);
double update_detsig_niso(double *lambda, double *psi, double log_detpsi, int p, int q);

// test for convergence (Aitken acceleration)
int converge_test(double *l, double *at, double *max_log_dens, double *log_dens, int N, int it, int G, double tol);

// matrix/vector operations
void generate_identity(int n, double *matrix);
void GaussJordan(int n, double *MATRIX, double *INVERSE, double *det);
void mx_mult(int m, int n, int q, double *a, double *b, double *r);
void mx_mult_diag1(int m, int n, double *a, double *b, double *r);
void vec_mx_mult(int n, int q, double *a, double *b, double *r);
void mx_trans(int m, int n, double *a, double *r);

// woodbury trick
void woodburyA(double *lambda, double psi, int p, int q, double *temp, double *cp, double *result);
double woodburyB(double *x, double psi, double *mu, int p, int q, double *lvec, double *cvec, double *cp);
void woodbury2A(double *lambda, double *psi, int p, int q, double *temp, double *cp, double *result);
double woodbury2B(double *x, double *psi, double *mu, int p, int q, double *lvec, double *cvec, double *cp);

// other
double maximum_array(double *array, int k);
void known_z(int *class, double *z, int N, int G);
void get_lambda(double *lam_vec, double **lambda, int G, int p, int q);
void store_lambda(double *lam_vec, double **lambda, int G, int p, int q);
