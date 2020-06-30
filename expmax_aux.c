#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "allocate.h"
#include "expmax.h"
#include "choldc.h"

float Quad_form_choldc(float *dx, float **sighalf, int Nvar);
float Determinant(float **A, int N);
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
int Sort(int *rank, float *score, int N);

/******************** Auxiliary routines ****************************/

int Sort(int *rank, float *score, int N)
{
  int i, j, k; float max=-1000;
  for(i=0; i<N; i++)rank[i]=-1;
  for(i=0; i<N; i++){
    k=-1; int *r=rank;
    for(j=0; j<N; j++){
      if((*r<0)&&((k<0)||(score[j]>max))){
	max=score[j]; k=j;
      }
      r++;
    }
    rank[k]=i;
  }
}

float Quad_form_choldc(float *dx, float **sighalf, int Nvar)
{
  double sum=0; int j;
  Forward_substitution_f(dx, sighalf, dx, Nvar);
  for(j=0; j<Nvar; j++)sum+=dx[j]*dx[j];
  return(sum);
}

float Determinant(float **A, int N)
{
  if(N==2){
    return(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
  }else if(N==3){
    double det=
      A[0][0]*A[1][1]*A[2][2]-
      A[0][0]*A[1][2]*A[2][1]+
      A[0][1]*A[1][2]*A[2][0]-
      A[0][1]*A[1][0]*A[2][2]+
      A[0][2]*A[1][0]*A[2][1]-
      A[0][2]*A[1][1]*A[2][0];
    return(det);
  }else if(N>3){
    int i,j,j1,j2,i2;
    double det = 0, sign=1.0;
    float **m = Allocate_mat2_f(N-1, N-1);
    for (j1=0; j1<N; j1++) {
      for (i=1; i<N; i++) {
	i2=i-1; j2=0;
	for (j=0; j<N; j++) {
	  if (j == j1)continue;
	  m[i2][j2] = A[i][j];
	  j2++;
	}
      }
      det +=  sign * A[0][j1] * Determinant(m,N-1);
      sign = -sign;
    }
    Empty_matrix_f(m, N-1);
    //printf("det = %.7f\n", det);
    return(det);
  }else if(N==1){
    return(A[0][0]);
  } else {
    printf("ERROR in determinant\n");
    return(0);
  }
}

double d_Determinant(double **A, int N)
{
  if(N==2){
    return(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
  }else if(N==3){
    double det=
      A[0][0]*A[1][1]*A[2][2]-
      A[0][0]*A[1][2]*A[2][1]+
      A[0][1]*A[1][2]*A[2][0]-
      A[0][1]*A[1][0]*A[2][2]+
      A[0][2]*A[1][0]*A[2][1]-
      A[0][2]*A[1][1]*A[2][0];
    return(det);
  }else if(N>3){
    int i,j,j1,j2,i2;
    double det = 0, sign=1.0;
    double **m = Allocate_mat2_d(N-1, N-1);
    for (j1=0; j1<N; j1++) {
      for (i=1; i<N; i++) {
	i2=i-1; j2=0;
	for (j=0; j<N; j++) {
	  if (j == j1)continue;
	  m[i2][j2] = A[i][j];
	  j2++;
	}
      }
      det +=  sign * A[0][j1] * d_Determinant(m,N-1);
      sign = -sign;
    }
    Empty_matrix_d(m, N-1);
    return(det);
  }else if(N==1){
    return(A[0][0]);
  } else {
    printf("ERROR in determinant\n");
    return(0);
  }
}

 
float **Transpose_matrix(float **m, int n1, int n2){
  float **mt=malloc(n2*sizeof(float *)); int i, j;
  for(i=0; i<n2; i++){
    float *a=malloc(n1*sizeof(float)); mt[i]=a;
    for(j=0; j<n1; j++){*a=m[j][i]; a++;}
  }
  return(mt);
}

float Dist2(float *x1, float *y1, int n)
{
  double d2=0; int j;
  float *x=x1, *y=y1, d;
  for(j=0; j<n; j++){d=*x-*y; d2+=d*d; x++; y++;}
  return(d2);
}

float Scalar_prod(float *u1, float *v1, int n){
  int i; float sum=0, *u=u1, *v=v1;
  for(i=0; i<n; i++){sum+=(*u)*(*v); u++; v++;}
  return(sum);
}

void Empty_variance_parameters(float *logdet, float ***sighalf,
			       float ***sigvec,
			       float  **sigval, float  **logval,
			       int *singular, int *print_sing,
			       float *logtau, int Nvar, int ncl)
{
  int j;
  for(j=0; j<ncl; j++){
    Empty_matrix_f(sighalf[j], Nvar);
    Empty_matrix_f(sigvec[j], Nvar);
    free(sigval[j]);
    free(logval[j]);
  }
  free(sighalf);
  free(sigvec);
  free(sigval);
  free(logval);
  free(singular);
  free(print_sing);
  free(logdet);
  free(logtau);
}

void Empty_Gaussian_parameters(float *tau, float **mu, float ***sig,
			       int Nvar, int ncl)
{
  int j;
  Empty_matrix_f(mu, ncl);
  for(j=0; j<ncl; j++)Empty_matrix_f(sig[j], Nvar);
  free(sig);
  free(tau);
}

/****************** Important functions! *******************/

int Process_parameters(float *logdet, float ***sighalf,  // Output
		       float ***sigvec,
		       float  **sigval, float  **logval,
		       int *singular, int *print_sing,
		       float *logtau, 
		       float ***sig, float *tau, // Input
		       int Nvar, int ncl, int iter)
{
  int i, j, k;
  for(k=0; k<ncl; k++){
    // Gaussian mixture
    float det=Determinant(sig[k], Nvar);
    singular[k]=choldc_f(sighalf[k], sig[k], Nvar);
    if((det<=0.0000001)||(singular[k]<0)){
      if(print_sing[k]==0){
	printf("WARNING ");
	if(singular[k]<0)printf("Choldc failed. ");
	printf("Singular matrix, cluster %d iter= %d Determinant= %.4f\n",
	       k, iter, det);
	print_sing[k]=1;
      }
      singular[k]=-1;
      f_Diagonalize(Nvar, sig[k], sigval[k], sigvec[k]);
      for(j=0; j<Nvar; j++){
	if(sigval[k][j]<lambda0)sigval[k][j]=lambda0;
	logval[k][j]=log(sigval[k][j]);
	sigval[k][j]=1./sqrt(sigval[k][j]);
      }
    }else{
      logdet[k]=log(det);
    }
    logtau[k]=log(tau[k]);
  }
  return(0);
}

double Compute_logGauss(float *dx, float *x,
			float *mu, float logdet,
			float **sighalf, float **sigvec,
			float *sigval, float *logval,
			int singular, int Nvar)
{
  double ee=0; int j;
  for(j=0; j<Nvar; j++)dx[j]=x[j]-mu[j];
  if(singular<0){
    for(j=0; j<Nvar; j++){
      float v=Scalar_prod(dx, sigvec[j], Nvar);
      v*=sigval[j];
      ee+=(v*v+logval[j]);
    }
  }else{
    Forward_substitution_f(dx, sighalf, dx, Nvar);
    ee=Scalar_prod(dx, dx, Nvar)+logdet;
  }
  return(-0.5*ee);
}

double Assign_clusters(int *cluster, float *Pstate, int *numcl, // Output 
		       float **log_Pxk, int N, int Nvar, int ncl)
{
  double lik=-0.5*(float)(N*Nvar)*log(2*3.14159);

  int i, k;
  for(k=0; k<ncl; k++)numcl[k]=0;
  for(i=0; i<N; i++){
    float *p=log_Pxk[i], pmax=p[0];
    int kmax=0;
    for(k=1; k<ncl; k++){
      if(p[k] > pmax){pmax=p[k]; kmax=k;}
    }
    Pstate[i]=exp(pmax);
    cluster[i]=kmax;
    numcl[kmax]++;
    lik+=pmax;
  }
  return(lik);
}

int Compactify_clusters(int *cluster, int *numcl, double *lik, // Output 
			float **Pxk, float **mu, float **x_SamVar,
			int N, int Nvar, int ncl)
{
  // Check direct interaction
  int nswitch=0, i, k, k1;
  float *dist2=malloc(ncl*sizeof(float));
  float **dc2=Allocate_mat2_f(ncl, ncl);
  for(k=0; k<ncl; k++){ 
    for(k1=0; k1<k; k1++){
      dc2[k][k1]=Dist2(mu[k], mu[k1], Nvar);
      dc2[k1][k]=dc2[k][k1];
    }
  }
  for(i=0; i<N; i++){
    int kmax=cluster[i], k2=-1;
    float *x=x_SamVar[i], dmin=-1;
    dist2[kmax]=Dist2(mu[kmax], x, Nvar);
    for(k=0; k<ncl; k++){
      if(k==kmax)continue;
      dist2[k]=Dist2(mu[k], x, Nvar);
      if((dist2[k]<dist2[kmax])&&(dc2[k][kmax]<dist2[kmax])){
	if((dmin<0)||(dist2[k]<dmin)){
	  dmin=dist2[k]; k2=k;
	}
      }
    }
    if(k2>=0){
      numcl[kmax]--; numcl[k2]++;
      float *p=Pxk[i];
      *lik-=log(p[kmax]);
      *lik+=log(p[k2]);
      cluster[i]=k2; nswitch++;
    }
  }
  printf("%d elements over %d changed cluster\n", nswitch, N);
  free(dist2); Empty_matrix_f(dc2, ncl);
}
