// Parameters internally set in expmax.c and also used in expmax_order.c :
extern float lambda0;   // Correction to singular matrix
extern int IT_MAX;      // Maximum number of iterations
extern float EPS;       // Threshold for convergence
extern int RMAX;        // Maximum number of initial conditions

// Main routines
float Expectation_maximization(int *cluster, float *Pstate, int ncl,
			       float *tau3, float **mu3,
			       float ***sig3,
			       float **x, int N, int Nvar, int first);
float Expectation_maximization_order(int *cluster, float *Pstate, int ncl,
				     float *tau3, float **corr3,
				     float **mu3, float ***sig3,
				     float **x, int N, int Nvar, int first);
// Common auxiliary routines
void Cluster_size(int *clus, int N, int ncl);
int Gaussian_parameters_1(float *tau, float **mu, float ***sig,
			  float **wki, float **x_SamVar,
			  int N, int Nvar, int ncl);
int Gaussian_parameters(float *tau, float **corr,
			float **mu, float ***sig,
			float **wik, float ***wikk1,
			float **x_SamVar,
			int N, int Nvar, int ncl, int iter);
int Initialize_parameters(float *tau, float **corr,
			  float **mu, float ***sig,
			  float **wki, float ***w2_ikk1,
			  float **x_SamVar, float *detsig,
			  int N, int Nvar, int ncl, int Round, int *cluster);

double Assign_clusters(int *cluster, float *Pstate, int *numcl, // Output
		       float **Pxk, int N, int Nvar, int ncl);
int Process_parameters(float *logdet, float ***sighalf,  // Output
		       float ***sigvec, float  **sigval, 
		       float  **logval,
		       int *singular, int *print_sing,
		       float *logtau, 
		       float ***sig, float *tau, // Input
		       int Nvar, int ncl, int iter);
double Compute_logGauss(float *dx, float *x,
			float *mu, float logdet,
			float **sighalf, float **sigvec,
			float *sigval, float *logval,
			int singular, int Nvar);
int Compactify_clusters(int *cluster, int *numcl, double *lik, // Output 
			float **w, float **mu, float **x_SamVar,
			int N, int Nvar, int ncl);

int Sort(int *rank, float *score, int N);
float Quad_form_choldc(float *dx, float **sighalf, int Nvar);
float Determinant(float **A, int N);
double d_Determinant(double **A, int N);
float Scalar_prod(float *u, float *v, int n);
void Empty_variance_parameters(float *logdet, float ***sighalf,
			       float ***sigvec,
			       float  **sigval, float  **logval,
			       int *singular, int *print_sing,
			       float *logtau, int Nvar, int ncl);
void Empty_Gaussian_parameters(float *tau, float **mu, float ***sig,
			       int Nvar, int ncl);

float **Transpose_matrix(float **m, int n1, int n2);
float Dist2(float *x1, float *y1, int n);
extern int COMPACT;
extern int ORDER;
