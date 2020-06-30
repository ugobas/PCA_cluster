#include <stdlib.h>
#include <math.h>
#include "allocate.h"
#include "expmax.h"
#include <stdio.h>

// Parameters of expectation maximization:
int RMAX=8;          // Maximum number of initial conditions
int IT_MAX=200;      // Maximum number of iterations
float EPS=0.000002;  // Threshold for convergence
float lambda0=0.0001; // Correction to singular matrix 

static int Gaussian_parameters(float *tau, float **mu, float ***sig,
			       float **wki, float **x_SamVar,
			       int N, int Nvar, int ncl);
static int Compute_weights(float **wki, float **log_Pxk, 
			   float *tau, float **mu,
			   float ***sig, float *logdet, float ***sighalf,
			   float ***sigvec, float **sigval, float **logval,
			   int *singular, int *print_sing, float *logtau,
			   float **x_SamVar, int N, int Nvar, int ncl,
			   int iter);
static void Copy_parameters(float *tau1, float **mu1, float ***sig1,
			    float *tau2, float **mu2, float ***sig2,
			    int ncl, int Nvar);
static float Check_convergence(float *tau1, float **mu1, float ***sig1,
			       float *tau2, float **mu2, float ***sig2,
			       int ncl, int Nvar);
static int Print_parameters(int *numcl, float *tau1, float **mu1,
			    float ***sig1, int ncl, int Nvar, int end);

int *print_nan;

/****************** Main routine ************************/
float Expectation_maximization(int *cluster, int ncl,
			       // Variables for storing optimal parameters:
			       float *tau3, float **mu3,
			       float ***sig3,
			       float **x_VarSam,
			       int N, int Nvar)
{
  double lik, lik_max=-10000000000;
  int i, j, j1, k, iter;
  int Round=1, ROUND_MAX=Nvar;
  if(ROUND_MAX > RMAX)ROUND_MAX=RMAX;

  float **x_SamVar=Transpose_matrix(x_VarSam, Nvar, N);
  
  /**************** Allocate memory ************************/
  // Variables used for weight computation
  float *tau1=malloc(ncl*sizeof(float));
  float **mu1=Allocate_mat2_f(ncl, Nvar);
  float ***sig1=malloc(ncl*sizeof(float **));
  for(k=0; k<ncl; k++)sig1[k]=Allocate_mat2_f(Nvar, Nvar);
  // Variables for updates
  float *tau2=malloc(ncl*sizeof(float));
  float **mu2=Allocate_mat2_f(ncl, Nvar);
  float ***sig2=malloc(ncl*sizeof(float **));
  for(k=0; k<ncl; k++)sig2[k]=Allocate_mat2_f(Nvar, Nvar);
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
  float **wki=Allocate_mat2_f(ncl, N);
  float **log_Pxk=Allocate_mat2_f(N, ncl);
  int *clus_max=malloc(N*sizeof(int));
  // Variables for storing optimal results
  // float *tau3=malloc(ncl*sizeof(float));
  // float **mu3=Allocate_mat2_f(ncl, Nvar);
  // float ***sig3=malloc(ncl*sizeof(float **));
  // for(k=0; k<ncl; k++)sig3[k]=Allocate_mat2_f(Nvar, Nvar);

 start_split:
  // Iterations start here
  for(k=0; k<ncl; k++){
    print_sing[k]=0; print_nan[k]=0;
  }
  for(i=0; i<N; i++)wki[0][i]=1;
  Gaussian_parameters(tau1, mu1, sig1, wki, x_SamVar, N, Nvar, 1);
  if(ncl==1){
    Process_parameters(logdet, sighalf, sigvec, sigval, logval,
		       singular, print_sing, logtau,
		       sig1, tau1, Nvar, ncl, iter);
    Compute_weights(wki, log_Pxk, tau1, mu1, sig1, logdet, sighalf,
		    sigvec, sigval, logval, singular, print_sing, logtau,
		    x_SamVar, N, Nvar, ncl, iter);
    goto Compute_lik;
  }

  /*************** Initial parameters ***************/
  printf("STARTING EXPECTATION MAXIMIZATION\n");
  // Check if the covariance is singular and regularize if not
  float det=Determinant(sig1[0], Nvar);
  if(det<=0.0000001){
    printf("WARNING, regularizing singular initial covariance matrix ");
    f_Diagonalize(Nvar, sig1[0], sigval[0], sigvec[0]);
    int n=0;
    for(k=0; k<Nvar; k++){
      if(sigval[0][k]<lambda0){
	printf("Singular component %d, sigma^2=%.5f\n", k, sigval[0][k]);
	sigval[0][k]=lambda0; n++;
      }
    }
    printf("Det= %.6f  %d singular components < %.5f\n", det, n, lambda0);
    for(j=0; j<Nvar; j++){
      for(j1=0; j1<=j; j1++){
	double sum=0;
	for(k=0; k<Nvar; k++)
	  sum+=sigval[0][k]*sigvec[0][k][i]*sigvec[0][k][j];
	sig1[0][j][j1]=sum; sig1[0][j1][j]=sum; 
      }
    }
  }
  // Parameters are copied from global statistics
  for(k=1; k<ncl; k++){
    for(j=0; j<Nvar; j++){
      mu1[k][j]=mu1[0][j];
      for(j1=0; j1<Nvar; j1++){
	sig1[k][j][j1]=sig1[0][j][j1];
	sig1[k][j1][j]=sig1[0][j][j1];
      }
    }
  }
  for(k=0; k<ncl; k++)tau1[k]=1./(float)ncl;
  // variable j=1...Round is distributed in [mu-2.5sigma,mu+2.5sigma]
  //for(j=0; j<Round; j++){
  j=Round-1;
  int l=ncl/2, l0=ncl-2*l, lk=l; k=0;
  float d=0.5*sig1[0][j][j]/(float)(l*l), mu=mu1[0][j];
  while(k<l){mu1[k][j]=mu-lk*lk*d; k++; lk--;}
  if(l0){mu1[k][j]=mu; k++;} lk=1;
  while(k<ncl){mu1[k][j]=mu+lk*lk*d; k++; lk++;}
  //}
  // Print initial parameters
  Print_parameters(numcl, tau1, mu1, sig1, ncl, Nvar, 0);
  /******************* End initial parameters ******************/

  /**** Expectation maximization algorithm for splitting N clusters ****/
  float err; int converged=1;
  for(iter=0; iter<IT_MAX; iter++){
    Process_parameters(logdet, sighalf, sigvec, sigval, logval,
		       singular, print_sing, logtau,
		       sig1, tau1, Nvar, ncl, iter);
    Compute_weights(wki, log_Pxk, tau1, mu1, sig1, logdet, sighalf,
		    sigvec, sigval, logval, singular, print_sing, logtau,
		    x_SamVar, N, Nvar, ncl, iter);
    //printf("%.0f\n", Assign_clusters(cluster, numcl, Pxk, N, Nvar, ncl));
    Gaussian_parameters(tau2, mu2, sig2, wki, x_SamVar, N, Nvar, ncl);
    err=Check_convergence(tau1, mu1, sig1, tau2, mu2, sig2, ncl, Nvar);
    Copy_parameters(tau1, mu1, sig1, tau2, mu2, sig2, ncl, Nvar);
    if(err<EPS)break;
  }
  if(iter==IT_MAX){
    printf("WARNING, Expectation maximization failed to converge,");
    converged=0;
  }else{
   printf("Expectation maximization converged,");
  }
  printf(" Attempt=%d iter= %d", Round, iter);
  printf(" error= %.8f Threshold=%.8f\n", err, EPS);

 Compute_lik:
  // Compute clusters
  lik=Assign_clusters(cluster, numcl, log_Pxk, N, Nvar, ncl);
  Print_parameters(numcl, tau1, mu1, sig1, ncl, Nvar, 1);
  printf("likelihood= %.0f\n", lik);

  if((Round==1)||(lik > lik_max)){
    lik_max=lik;
    for(i=0; i<N; i++)clus_max[i]=cluster[i];
    Copy_parameters(tau3, mu3, sig3, tau1, mu1, sig1, ncl, Nvar); 
  }
  Round++;
  if((Round <= ROUND_MAX)&&(ncl > 1)){   //(converged==0)&&
    goto start_split;
  }

  for(i=0; i<N; i++)cluster[i]=clus_max[i];
  Print_parameters(numcl, tau3, mu3, sig3, ncl, Nvar, 1);
  printf("likelihood= %.0f\n", lik_max);

  if(COMPACT){
    Compactify_clusters(cluster,numcl,&lik,log_Pxk,mu3,x_SamVar,N,Nvar,ncl);
  }

  // Empty
  Empty_matrix_f(x_SamVar, N);
  Empty_matrix_f(wki, ncl);
  Empty_matrix_f(log_Pxk, N);
  free(numcl);
  free(clus_max);
  free(print_nan);
  Empty_Gaussian_parameters(tau1, mu1, sig1, Nvar, ncl);
  Empty_Gaussian_parameters(tau2, mu2, sig2, Nvar, ncl);
  Empty_variance_parameters(logdet, sighalf, sigvec, sigval, logval,
			    singular, print_sing, logtau, Nvar, ncl);
  return(lik);
}

static int Gaussian_parameters(float *tau, float **mu, float ***sig,
			       float **wki, float **x_SamVar,
			       int N, int Nvar, int ncl)
{
  int i, j, j1, k; 
  double *x1=malloc(Nvar*sizeof(double));
  double **x2=Allocate_mat2_d(Nvar, Nvar);

  for(k=0; k<ncl; k++){
    for(j=0; j<Nvar; j++){
      x1[j]=0; for(j1=0; j1<=j; j1++)x2[j][j1]=0;
    }
    // Sum over all samples
    double sum=0; int nnan=0;
    float *w=wki[k];
    for(i=0; i<N; i++){
      if(isnan(*w)){*w=0; nnan++;}
      float *x=x_SamVar[i], wx;
      for(j=0; j<Nvar; j++){
	wx=(*w)*x[j]; x1[j]+=wx;
	for(j1=0; j1<=j; j1++)x2[j][j1]+=wx*x[j1];
      }
      sum+=(*w); w++;
    }
    // Warn if not a number 
    if((nnan)&&(print_nan[k]==0)){
      printf("WARNING, cl. %d %d elements are nan\n", k, nnan);
      print_nan[k]=1;
    }
    if((sum<=0)||(isnan(sum))){
      printf("WARNING, cl. %d sum of weights= %f nan in w: ", k, sum);
      printf("\n");
      return(-1);
    }
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

  return(0);
}

static int Compute_weights(float **wki, float **log_Pxk,
			   float *tau, float **mu, float ***sig,
			   float *logdet, float ***sighalf, float ***sigvec,
			   float **sigval, float **logval,
			   int *singular, int *print_sing, float *logtau,
			   float **x_SamVar, int N, int Nvar, int ncl, int iter)
{
  Process_parameters(logdet, sighalf, sigvec, sigval, logval,
		     singular, print_sing, logtau,
		     sig, tau, Nvar, ncl, iter);

  int i, k, failed=0;
  float *dx=malloc(Nvar*sizeof(float));
  for(k=0; k<ncl; k++){
    int nnan=0, sing=singular[k];
    float *m=mu[k], ld=logdet[k], **sh=sighalf[k], t=tau[k], logw;
    float **svec=sigvec[k], *sval=sigval[k], *lval=logval[k], lt=logtau[k];
    for(i=0; i<N; i++){
      logw=Compute_logGauss(dx,x_SamVar[i],m,ld,sh,svec,sval,lval,sing,Nvar);
      if(isnan(logw)){
	nnan++; logw=-100;
      }
      log_Pxk[i][k]=logw+lt;
    }
    if(nnan){
      printf("WARNING, cluster %d has %d weights nan\n", k, nnan);
      printf("tau: %.2f  logdet: %.2f Determinant= %.4f\n",
	     t, ld, Determinant(sig[k], Nvar));
      if(nnan>10)failed=1;
    }
  }
  float *p=malloc(ncl*sizeof(float));
  for(i=0; i<N; i++){
    double sum=0;
    float *lp=log_Pxk[i];
    for(k=0; k<ncl; k++){p[k]=exp(lp[k]); sum+=p[k];}
    for(k=0; k<ncl; k++)wki[k][i]=p[k]/sum; 
  }
  free(dx); free(p);

  return(0);
}

static void Copy_parameters(float *tau1, float **mu1, float ***sig1,
			    float *tau2, float **mu2, float ***sig2,
			    int ncl, int Nvar)
{
  int k, j, j1;
  for(k=0; k<ncl; k++){
    tau1[k]=tau2[k];
    float *m1=mu1[k], *m2=mu2[k], **s1=sig1[k], **s2=sig2[k];
    for(j=0; j<Nvar; j++){
      m1[j]=m2[j];
      for(j1=0; j1<Nvar; j1++)s1[j][j1]=s2[j][j1];
    }
  }
}

static float Check_convergence(float *tau1, float **mu1, float ***sig1,
			       float *tau2, float **mu2, float ***sig2,
			       int ncl, int Nvar)
{
  double sum=0; int i, j, j1, k;
  for(k=0; k<ncl; k++){
    sum += fabs(tau1[k]-tau2[k]);
    float *m1=mu1[k], *m2=mu2[k], **s1=sig1[k], **s2=sig2[k];
    for(j=0; j<Nvar; j++){
      sum += fabs(m1[j]-m2[j])/fabs(m1[j]+m2[j]);
      for(j1=j; j1<Nvar; j1++){
	sum += fabs(s1[j][j1]-s2[j][j1])/fabs(s1[j][j1]+s2[j][j1]);
      }
    }
  }
  int Npar = ncl*(1+Nvar+Nvar*(Nvar+1)/2)-1;
  return(sum/Npar);
}

static int Print_parameters(int *numcl, float *tau1,
			    float **mu1, float ***sig1,
			    int ncl, int Nvar, int end)
{
  int k, j, i;
  for(k=0; k<ncl; k++){
    if(end){
      printf("Clus%d %4d samples, ", k+1, numcl[k]);
    }else{
      printf("Clus%d Initial para: ", k+1);
    }
    printf("tau= %.2f mu= ", tau1[k]);
    for(j=0; j<Nvar; j++)printf("%6.3f ", mu1[k][j]);
    printf(" sig= ");
    for(j=0; j<Nvar; j++)printf("%5.3f ", sig1[k][j][j]);
    /*printf(" ");
    for(j=0; j<Nvar; j++)
    for(i=0; i<j; i++)printf("%6.3f ", sig1[k][j][i]);*/
    printf("\n");
  }
}


