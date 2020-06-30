#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "allocate.h"
#include "expmax.h"
#include "average_linkage.h"

// Internal parameters defined in expmax.h and assigned in expmax.c:
// int RMAX=5;          // Maximum number of initial conditions
// int IT_MAX=200;      // Maximum number of iterations
// float EPS=0.000002;  // Threshold for convergence
// float lambda0=0.0001; // Correction to singular matrix 


int Stationary(float *P_i, float **P_ij, int N);
static double Viterbi(int *state, int *num,
		      float **lgik, float *tau, float **corr,
		      int N, int Nvar, int ncl);
static int
Forward_backward(float **w1_ik, float ***w2_ik, float **lg_ik, // Output
		 float **x_SamVar, float *tau1, float **corr1,
		 float **mu1, float ***sig1, int N, int Nvar, int ncl, // input
		 float **wf_ik, float **g_ik, float *logtau, float *logdet,
		 float ***sighalf, float ***sigvec, float **sigval,
		 float **logval, int *singular, int *print_sing, int iter);

int Gaussian_parameters(float *tau, float **corr,
			float **mu, float ***sig,
			float **wik, float ***wikk1,
			float **x_SamVar,
			int N, int Nvar, int ncl, int iter);
static int Compute_weights(float **wik, float ***wikk1, float **lgik,
			   float *tau, float **corr, float **mu, float ***sig,
			   float *logdet, float ***sighalf, float ***sigvec,
			   float **sigval, float **logval,
			   int *singular, int *print_sing, float *logtau,
			   float **x_SamVar,int N, int Nvar, int ncl, int iter);
static void Copy_parameters(float *tau1, float **corr1,
			    float **mu1, float ***sig1,
			    float *tau2, float **corr2,
			    float **mu2, float ***sig2,
			    int ncl, int Nvar);
static float Check_convergence(float *tau1, float **corr1,
			       float **mu1, float ***sig1,
			       float *tau2, float **corr2,
			       float **mu2, float ***sig2,
			       int ncl, int Nvar);
static int Print_parameters(int *numcl, float *tau, float **corr,
			    float **mu, float ***sig, int ncl,
			    int Nvar, int end);
void Copy_mu(float **Mu, float **mu1, int ncl, int Nvar);
float Compare_mu(float ***Mu, int nround,  int ncl, int Nvar);
float Cosine(float *v1, float *v2, int N);

int *print_nan;

/****************** Main routine ************************/
float Expectation_maximization_order(int *clustout, float *Pstate_out, int ncl,
				     float *tau3, float **corr3,
				     float **mu3, float ***sig3,
				     float **x_VarSam,
				     int N, int Nvar, int first)
{
  double lik, lik_max=-10000000000;
  int i, j, j1, k, k1, iter, Round=1;
  int ROUND_MAX=RMAX;

  float **x_SamVar=Transpose_matrix(x_VarSam, Nvar, N);
  
  /**************** Allocate memory ************************/
  // Variables used for weight computation
  float *tau1=malloc(ncl*sizeof(float));
  float **corr1=Allocate_mat2_f(ncl, ncl);
  float **mu1=Allocate_mat2_f(ncl, Nvar);
  float ***sig1=malloc(ncl*sizeof(float **));
  for(k=0; k<ncl; k++)sig1[k]=Allocate_mat2_f(Nvar, Nvar);
  // Variables for updates
  float *tau2=malloc(ncl*sizeof(float));
  float **corr2=Allocate_mat2_f(ncl, ncl);
  float **mu2=Allocate_mat2_f(ncl, Nvar);
  float ***sig2=malloc(ncl*sizeof(float **));
  for(k=0; k<ncl; k++)sig2[k]=Allocate_mat2_f(Nvar, Nvar);
  // Variables for storing optimal results are passed from calling program

  // Variables for processing variance
  float ***sighalf=malloc(ncl*sizeof(float **));
  float *logdet=malloc(ncl*sizeof(float));
  float ***sigvec=malloc(ncl*sizeof(float **));
  float  **sigval=Allocate_mat2_f(ncl, Nvar);
  float  **logval=Allocate_mat2_f(ncl, Nvar);
  int *singular=malloc(ncl*sizeof(int));
  int *print_sing=malloc(ncl*sizeof(int));
  print_nan=malloc(ncl*sizeof(int));
  float *logtau=malloc(ncl*sizeof(float));
  for(k=0; k<ncl; k++){
    sighalf[k]=Allocate_mat2_f(Nvar, Nvar);
    sigvec[k]=Allocate_mat2_f(Nvar, Nvar);
  }
  // Variables for weights and clusters
  int *numcl=malloc(ncl*sizeof(int));
  float **w1_ik=Allocate_mat2_f(N, ncl);
  float ***w2_ikk1=malloc(N*sizeof(float **));
  for(i=0; i<N; i++)w2_ikk1[i]=Allocate_mat2_f(ncl, ncl);
  float **wf_ik=Allocate_mat2_f(N, ncl);
  float **lg_ik=Allocate_mat2_f(N, ncl);
  float **g_ik=Allocate_mat2_f(N, ncl);
  int *clus_max=malloc(N*sizeof(int));
  int *cluster=malloc(N*sizeof(int));

  float ***Mu=malloc(ROUND_MAX*sizeof(float **));
  for(k=0; k<ROUND_MAX; k++)Mu[k]=Allocate_mat2_f(ncl, Nvar);

  // Compute variance for ranking previous clusters
  float *detsig=malloc((ncl-1)*sizeof(float));
  for(k=0; k<(ncl-1); k++)detsig[k]=Determinant(sig3[k], Nvar);

  // ncl different initializations
  int Rmax=ncl-1;
  if(Rmax==0)Rmax=1;
  if(Rmax > ROUND_MAX)Rmax=ROUND_MAX;
  int *clustini=clustout; if(first)clustini=NULL;

  for(Round=0; Round<Rmax; Round++){

    /****************** Iterations start here ****************/
    for(k=0; k<ncl; k++){
      print_sing[k]=0; print_nan[k]=0;
    }

    /*************** Initial parameters ***************/
    Initialize_parameters(tau1, corr1, mu1, sig1, w1_ik, w2_ikk1, x_SamVar,
			  detsig, N, Nvar, ncl, Round, clustini);
    Print_parameters(numcl, tau1, corr1, mu1, sig1, ncl, Nvar, 0);
    if(ncl==1)goto Compute_lik;

    /**** Expectation maximization algorithm for splitting N clusters ****/
    float err; int converged=1;
    for(iter=0; iter<IT_MAX; iter++){
      Process_parameters(logdet, sighalf, sigvec, sigval, logval,
			 singular, print_sing, logtau,
			 sig1, tau1, Nvar, ncl, iter);
      Forward_backward(w1_ik, w2_ikk1, NULL,
		       x_SamVar, tau1, corr1, mu1, sig1, N, Nvar, ncl,
		       wf_ik, g_ik, logtau, logdet, sighalf, sigvec,
		       sigval, logval, singular, print_sing, iter);
      Gaussian_parameters(tau2, corr2, mu2, sig2,
			  w1_ik, w2_ikk1, x_SamVar, N, Nvar, ncl, iter);
      err=Check_convergence(tau1, corr1, mu1, sig1,
			    tau2, corr2, mu2, sig2, ncl, Nvar);
      //Print_parameters(numcl, tau2, corr2, mu2, sig2, ncl, Nvar, 0);
      Copy_parameters(tau1,corr1,mu1,sig1, tau2,corr2,mu2,sig2, ncl, Nvar);
      if(err<EPS)break;
    }
    if(iter==IT_MAX){
      printf("WARNING, HMM failed to converge,");
      converged=0;
    }else{
      printf("HMM converged,");
    }
    printf(" Round=%d iter= %d", Round, iter);
    printf(" error= %.8f Threshold=%.8f\n", err, EPS);

  Compute_lik:
    // Compute clusters
    Forward_backward(w1_ik, w2_ikk1, lg_ik,
		     x_SamVar, tau1, corr1, mu1, sig1, N, Nvar, ncl,
		     wf_ik, g_ik, logtau, logdet, sighalf, sigvec,
		     sigval, logval, singular, print_sing, iter);
    lik=Viterbi(cluster, numcl, lg_ik, logtau, corr1, N, Nvar, ncl);
    Print_parameters(numcl, tau1, corr1, mu1, sig1, ncl, Nvar, 1);
    Copy_mu(Mu[Round], mu1, ncl, Nvar);
    printf("likelihood/N= %.3f\n", lik/N);
    if((Round==0)||(lik > lik_max)){
      lik_max=lik;
      for(i=0; i<N; i++){
	clus_max[i]=cluster[i];
	Pstate_out[i]=w1_ik[i][cluster[i]];
      }
      Copy_parameters(tau3,corr3,mu3,sig3,tau1,corr1,mu1,sig1,ncl,Nvar); 
    }
  }
  free(cluster);
  free(detsig);

  for(i=0; i<N; i++)clustout[i]=clus_max[i];
  printf("All rounds: likelihood/N= %.3f\n", lik_max/N);
  Print_parameters(numcl, tau3, corr3, mu3, sig3, ncl, Nvar, 1);

  /*if(COMPACT){
    Compactify_clusters(cluster, numcl, &lik, Pxk, mu3, x_SamVar, N, Nvar, ncl);
    }*/
  if(Rmax>1){
    double sim=Compare_mu(Mu, Round, ncl, Nvar);
    printf("Similarity between %d different initial conditions: %.3f\n",
	   Rmax, sim);
  }

  // Empty
  Empty_matrix_f(x_SamVar, N);
  Empty_matrix_f(wf_ik, N);
  Empty_matrix_f(lg_ik, N);
  Empty_matrix_f(g_ik, N);
  Empty_matrix_f(w1_ik, N);
  for(i=0; i<N; i++)Empty_matrix_f(w2_ikk1[i], ncl);
  free(w2_ikk1);
  free(numcl);
  free(clus_max);
  free(print_nan);
  Empty_Gaussian_parameters(tau1, mu1, sig1, Nvar, ncl);
  Empty_Gaussian_parameters(tau2, mu2, sig2, Nvar, ncl);
  Empty_variance_parameters(logdet, sighalf, sigvec, sigval, logval,
			    singular, print_sing, logtau, Nvar, ncl);
  Empty_matrix_f(corr1, ncl);
  Empty_matrix_f(corr2, ncl);
  for(k=0; k<ROUND_MAX; k++)Empty_matrix_f(Mu[k], ncl);
  free(Mu);
  //Empty_Gaussian_parameters(tau3,mu3,sig3, Nvar, ncl);
  //Empty_matrix_f(corr3, ncl);
  ///////////////////////////////////////////////////

  return(lik_max);
}

int Gaussian_parameters(float *tau, float **corr,
			float **mu, float ***sig,
			float **wik, float ***wikk1,
			float **x_SamVar,
			int N, int Nvar, int ncl, int iter)
{
  int i, j, j1, k, empty=0; 
  double *x1=malloc(Nvar*sizeof(double));
  double **x2=Allocate_mat2_d(Nvar, Nvar);

  for(k=0; k<ncl; k++){
    for(j=0; j<Nvar; j++){
      x1[j]=0; for(j1=0; j1<=j; j1++)x2[j][j1]=0;
    }
    // Sum over all samples
    double sum=0; int nnan=0;
    for(i=0; i<N; i++){
      float w=wik[i][k];
      float *x=x_SamVar[i], wx;
      if(isnan(w)){w=0; nnan++;}
      for(j=0; j<Nvar; j++){
	wx=w*x[j]; x1[j]+=wx;
	for(j1=0; j1<=j; j1++)x2[j][j1]+=wx*x[j1];
      }
      sum+=w;
    }
    // Warn if not a number 
    if((nnan)&&(print_nan[k]==0)){
      printf("WARNING, cl. %d %d elements are nan iter=%d\n", k, nnan, iter);
      print_nan[k]=1;
    }
    if((sum<=0)||(isnan(sum))){
      printf("WARNING, cl. %d sum of weights= %f, iter=%d\n", k, sum, iter);
      empty=1;
    }
    if(empty)printf("cl. %d sum of weights= %f, iter=%d\n", k, sum, iter);
    // Compute averages
    tau[k]=sum/N;
    float **s=sig[k];
    for(j=0; j<Nvar; j++){
      x1[j]/=sum;  mu[k][j]=x1[j];
      for(j1=0; j1<=j; j1++){
	s[j][j1]=x2[j][j1]/sum-x1[j]*x1[j1];
	if(j!=j1)s[j1][j]=s[j][j1];
      }
    }
  }
  free(x1); Empty_matrix_d(x2, Nvar);

  if(corr){

    // Compute correlations
    int k1;
    double **cc=Allocate_mat2_d(ncl, ncl), sum;
    for(i=0; i<N-1; i++){
      for(k1=0; k1<ncl; k1++){
	double *ck=cc[k1]; float *wk=wikk1[i][k1];
	for(k=0; k<ncl; k++)ck[k]+=wk[k];
      }
    }
    for(k1=0; k1<ncl; k1++){
      double *ck=cc[k1];
      sum=0; for(k=0; k<ncl; k++)sum+=ck[k];
      float *c2=corr[k1];
      for(k=0; k<ncl; k++)c2[k]=ck[k]/sum;
    }
    Empty_matrix_d(cc, ncl);

    // Compute tau as stationary distribution of corr
    Stationary(tau, corr, ncl);
  }
  return(0);
}

int Stationary(float *P_i, float **P_ij, int N)
{
  int ITMAX=100; float thr=0.00001*N;
  float *P2_i=malloc(N*sizeof(float));
  int it, i, j;
  for(i=0; i<N; i++)P_i[i]=1./(float)N;
  for(it=0; it<ITMAX; it++){
    double sum=0, err=0, p;
    for(i=0; i<N; i++){
      p=0; for(j=0; j<N; j++)p+=P_i[j]*P_ij[j][i];
      P2_i[i]=p; sum+=p;
    }
    for(i=0; i<N; i++){
      p=P2_i[i]/sum;
      err+=fabs(p-P_i[i]);
      P_i[i]=p;
    }
    if(err < thr)break;
  }
  //printf("Stationary distribution converged in %d steps, error= %.6f\n");
  //for(i=0; i<N; i++)printf("%.2f ", P_i[i]); printf("\n");
  free(P2_i);
}
      

static int
Forward_backward(float **w1_ik, float ***w2_ikk1, float **lg_ik, // Output
		 float **x_SamVar, float *tau, float **corr,
		 float **mu, float ***sig, int N, int Nvar, int ncl, // input
		 float **wf_ik, float **g_ik, float *logtau, float *logdet,
		 float ***sighalf, float ***sigvec, float **sigval,
		 float **logval, int *singular, int *print_sing, int iter)
{

  // Preparation
  float *dx=malloc(Nvar*sizeof(float));
  float *z_i=malloc(N*sizeof(float));
  int *nnan=malloc(ncl*sizeof(int));
  int i, i1, k, k1, failed=0;
  int zero=0, zerot=0;
  float epsg=0.000001;
  for(k=0; k<ncl; k++)nnan[k]=0;
  Process_parameters(logdet, sighalf, sigvec, sigval, logval,
		     singular, print_sing, logtau,
		     sig, tau, Nvar, ncl, iter);

  /********************  Forward step  ************************/
  float *nu=malloc(ncl*sizeof(float));
  for(k=0; k<ncl; k++)nu[k]=tau[k];
  for(i=0; i<N; i++){
    float *x=x_SamVar[i], *w=wf_ik[i];
    float g, *lg, *gi=g_ik[i], logw;
    if(lg_ik)lg=lg_ik[i];
    double sum=0;
    for(k=0; k<ncl; k++){
      logw=Compute_logGauss(dx, x, mu[k], logdet[k], sighalf[k], sigvec[k],
			    sigval[k], logval[k], singular[k], Nvar);

      if(isnan(logw)){logw=-100; nnan[k]++;}
      g=exp(logw);
      if(g< epsg)g=epsg;
      gi[k]=g; if(lg_ik)lg[k]=logw;
      g*=nu[k]; w[k]=g; sum+=g;
    }
    if(sum==0){
      if(zero==0){
	printf("First zero weight, i=%d zerot=%d sum=%.f\n", i,zerot, sum);
      }
      float p=1./ncl; for(k=0; k<ncl; k++)w[k]=p; sum=epsg;
      zero++;
    }
    z_i[i]=sum;
    for(k=0; k<ncl; k++)nu[k]=0;
    for(k=0; k<ncl; k++){
      w[k]/=sum; float *q=corr[k];
      for(k1=0; k1<ncl; k1++)nu[k1]+=w[k]*q[k1];
    }
  }
  free(dx);

  if(zerot){
    printf("WARNING, %d prefactors are zero, iter=%d\n", zerot, iter);
  }
  if(zero){
    printf("WARNING, %d weights are zero, iter=%d\n", zero, iter);
  }
  for(k=0; k<ncl; k++){
    if(nnan[k]){
      printf("WARNING, cluster %d has %d weights nan iter=%d\n",
	     k, nnan[k], iter);
      printf("tau: %.2f  logdet: %.2f Determinant= %.4f\n",
	     tau[k], logdet[k], Determinant(sig[k], Nvar));
      if(nnan[k]>10)failed=1;
    }
  }
  free(nnan);
  free(nu);

  /********************  Backward step  ************************/
  // Initialize
  i=N-1;
  float zz=1./z_i[i], *w=w1_ik[i], *wf=wf_ik[i];
  float *beta1=malloc(ncl*sizeof(float));
  float *beta2=malloc(ncl*sizeof(float)), b;
  double norm=0;
  for(k=0; k<ncl; k++){
    beta1[k]=zz; norm+=wf[k];
  }
  for(k=0; k<ncl; k++)w[k]=wf[k]/norm;
  // Iterate
  for(i=N-2; i>=0; i--){
    w=w1_ik[i]; wf=wf_ik[i];
    zz=1./z_i[i]; i1=i+1; norm=0;
    float *g1=g_ik[i1], **w2=w2_ikk1[i];
    for(k=0; k<ncl; k++){
      double sum=0; float *q=corr[k];
      for(k1=0; k1<ncl; k1++){
	float Trans=q[k1]*g1[k1]*beta1[k1];
	w2[k][k1]=wf[k]*Trans;
	sum+=Trans;
      }
      b=sum*zz;
      beta2[k]=b;
      w[k]=wf[k]*b;
      norm+=w[k];
    }
    for(k=0; k<ncl; k++){
      beta1[k]=beta2[k];
      w[k]/=norm;
    }
  }
  free(beta1); free(beta2);
  free(z_i);
}

static float Check_convergence(float *tau1, float **corr1,
			       float **mu1, float ***sig1,
			       float *tau2, float **corr2,
			       float **mu2, float ***sig2,
			       int ncl, int Nvar)
{
  double sum=0; int i, j, j1, k, k1;
  for(k=0; k<ncl; k++){
    for(k1=0; k1<ncl; k1++){
      sum += fabs(corr1[k][k1]-corr2[k][k1]);
    }
    float *m1=mu1[k], *m2=mu2[k], **s1=sig1[k], **s2=sig2[k];
    for(j=0; j<Nvar; j++){
      sum += fabs(m1[j]-m2[j])/fabs(m1[j]+m2[j]);
      sum += fabs(s1[j][j]-s2[j][j])/fabs(s1[j][j]+s2[j][j]);
      for(j1=0; j1<j; j1++){
	sum += fabs(s1[j][j1]-s2[j][j1])/fabs(s1[j][j1]+s2[j][j1]);
      }
    }
  }
  int Npar = ncl*((ncl-1)+(Nvar*(Nvar+1))/2+(Nvar+1));
  return(sum/Npar);
}

static int Print_parameters(int *numcl,
			    float *tau, float **corr,
			    float **mu, float ***sig,
			    int ncl, int Nvar, int end)
{
  int k, j, i;
  for(k=0; k<ncl; k++){
    if(end){
      printf("Clus%d %4d samples, ", k+1, numcl[k]);
    }else{
      printf("Clus%d Initial para: ", k+1);
    }
    printf("tau= %.2f mu= ", tau[k]);
    for(j=0; j<Nvar; j++)printf("%6.3f ", mu[k][j]);
    printf(" sig= ");
    for(j=0; j<Nvar; j++)printf("%5.3f ", sig[k][j][j]);
    printf(" corr= ");
    for(j=0; j<ncl; j++)printf("%4.2f ", corr[k][j]);
    printf("\n");
  }
}

static void Copy_parameters(float *tau1, float **corr1,
			    float **mu1, float ***sig1,
			    float *tau2, float **corr2,
			    float **mu2, float ***sig2,
			    int ncl, int Nvar)
{
  int k, j, j1;
  for(k=0; k<ncl; k++){
    tau1[k]=tau2[k];
    float *m1=mu1[k], *m2=mu2[k], **s1=sig1[k], **s2=sig2[k];
    float *c1=corr1[k], *c2=corr2[k];
    for(j=0; j<ncl; j++)c1[j]=c2[j];
    for(j=0; j<Nvar; j++){
      m1[j]=m2[j];
      for(j1=0; j1<Nvar; j1++)s1[j][j1]=s2[j][j1];
    }
  }
}

double Viterbi(int *state, int *num, float **lgik,
	       float *logtau, float **corr,
	       int N, int Nvar, int ncl)
{  
  double lik=-0.5*(float)(N*Nvar)*log(2*3.14159);
  int i, k, k1, k1max;
  for(k=0; k<ncl; k++)num[k]=0;

  // Forward step
  float **sik=Allocate_mat2_f(N, ncl), s, smax;
  float **lcorr=Allocate_mat2_f(ncl, ncl);
  for(k=0; k<ncl; k++){
    sik[0][k]=logtau[k]+lgik[0][k];
    for(k1=0; k1<ncl; k1++)lcorr[k][k1]=log(corr[k][k1]);
  }
  for(i=1; i<N; i++){
    float *si=sik[i-1];
    for(k=0; k<ncl; k++){
      for(k1=0; k1<ncl; k1++){
	s=si[k1]+lcorr[k1][k]+lgik[i][k];
	if((s>smax)||(k1==0)){smax=s; k1max=k1;}
      }
      sik[i][k]=smax;
    }
  }

  // Backward step
  int kmax, i1; float *si;
  i=N-1; si=sik[i];
  for(k=0; k<ncl; k++){
    if((si[k]>smax)||(k==0)){smax=si[k]; kmax=k;}
  }
  num[kmax]++;
  state[i]=kmax;
  for(i=N-2; i>=0; i--){
    si=sik[i];
    i1=i+1;
    k1=state[i1];
    for(k=0; k<ncl; k++){
      s=si[k]+lcorr[k][k1]+lgik[i1][k1];
      if((s>smax)||(k==0)){smax=s; kmax=k;}
    }
    num[kmax]++;
    state[i]=kmax;
    lik+=lcorr[kmax][k1]+lgik[i1][k1];
  }
  lik+=logtau[kmax]+lgik[0][kmax];

  Empty_matrix_f(sik, N);
  Empty_matrix_f(lcorr, ncl);
  return(lik);
}

static int Compute_weights(float **wik, float ***wikk1, float **lgik,
			   float *tau, float **corr, float **mu, float ***sig,
			   float *logdet, float ***sighalf, float ***sigvec,
			   float **sigval, float **logval,
			   int *singular, int *print_sing, float *logtau,
			   float **x_SamVar, int N, int Nvar, int ncl, int iter)
{
  Process_parameters(logdet, sighalf, sigvec, sigval, logval,
		     singular, print_sing, logtau,
		     sig, tau, Nvar, ncl, iter);

  int i, i1, k, k1, failed=0;
  float *dx=malloc(Nvar*sizeof(float));
  int *nnan=malloc(ncl*sizeof(int));
  for(k=0; k<ncl; k++)nnan[k]=0;

  for(i=0; i<N; i++){
    float *x=x_SamVar[i], *w=wik[i], *pk;
    float t, g, ww, *lg, logw;
    double sum=0;
    for(k=0; k<ncl; k++){
      logw=Compute_logGauss(dx, x, mu[k], logdet[k], sighalf[k], sigvec[k],
			    sigval[k], logval[k], singular[k], Nvar);
      g=exp(logw);
      if(isnan(g)){nnan[k]++; g=0;}
      if(lgik)lg[k]=logw+logtau[k];

      //Marginalize(tauc, corr, w, tau, k, i, ncl);
      float pik;
      if(i){
	float *w1=wik[i-1], **w2=wikk1[i-1]; pik=0;
	for(k1=0; k1<ncl; k1++){
	  float p=w1[k1]*corr[k1][k]*g;
	  w2[k1][k]=p;
	  pik+=p;
	}
      }else{
	pik=tau[k]*g;
      }
      w[k]=pik;
      sum+=pik;
    }
    if(sum>0){
      for(k=0; k<ncl; k++)w[k]/=sum;
    }else{
      float p=1./ncl;
      for(k=0; k<ncl; k++)w[k]=p;
    }
  }
  free(dx);

  for(k=0; k<ncl; k++){
    if(nnan[k]){
      printf("WARNING, cluster %d has %d weights nan iter=%d\n",
	     k, nnan[k], iter);
      printf("tau: %.2f  logdet: %.2f Determinant= %.4f\n",
	     tau[k], logdet[k], Determinant(sig[k], Nvar));
      if(nnan[k]>10)failed=1;
    }
  }
  free(nnan);

  return(0); // failed
}

void Copy_mu(float **Mu, float **mu1, int ncl, int Nvar){
  int k, j;
  for(k=0; k<ncl; k++){
    for(j=0; j<Nvar; j++)Mu[k][j]=mu1[k][j];
  }
}

float Compare_mu(float ***Mu, int nround,  int ncl, int Nvar)
{
  int i1, i2, k1, k2, j;
  float **Mu1, **Mu2;
  //int *coclus=malloc(ncl*sizeof(int));
  float *cosclus=malloc(ncl*sizeof(float)), c;
  double simsum=0; long norm=0, sum=0;
  /*float *w=malloc(ncl*sizeof(float));
  for(k1=0; k1<ncl; k1++)sum+=num[k1];
  for(k1=0; k1<ncl; k1++)w[k1]=num[k1]/sum;*/
  for(i1=0; i1<nround; i1++){
    Mu1=Mu[i1];
    for(i2=0; i2<i1; i2++){
      Mu2=Mu[i2];
      float sim=0; 
      for(k1=0; k1<ncl; k1++){
	cosclus[k1]=-2;
	for(k2=0; k2<ncl; k2++){
	  c=Cosine(Mu1[k1], Mu2[k2], Nvar);
	  if(c>cosclus[k1])cosclus[k1]=c;
	}
	sim+=cosclus[k1];
      }
      simsum+=sim/ncl; norm++;
    }
  }
  free(cosclus); //free(w);
  return(simsum/norm);
}
