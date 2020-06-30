#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "expmax.h"
#include "allocate.h"
#include "random3.h"
#include <time.h>

int *sort=NULL;
float **sigvec=NULL;
float *sigval=NULL;
RANDOMTYPE iran;

int *Rank_clus(float *detsig, int nc1);
void Compute_trans(float **trans, float ***w2_ikk1,
		   int *cluster, int ncl, int clus, int N);
unsigned long randomgenerator(void);

int Initialize_parameters(float *tau, float **trans,
			  float **mu, float ***sig,
			  float **w1_ik, float ***w2_ikk1,
			  float **x_SamVar, float *detsig,
			  int N, int Nvar, int ncl, int Round, int *cluster)
{
  int nc1=ncl-1, i, j, j1, k, k1, clus=0;

  if(ncl==1){
    for(i=0; i<N; i++)w1_ik[i][0]=1;
    if(w2_ikk1)for(i=0; i<N; i++)w2_ikk1[i][0][0]=1;
    if(trans)trans[0][0]=1;
    if(sigvec)Empty_matrix_f(sigvec, Nvar);
    if(sigval)free(sigval);
    sigvec=Allocate_mat2_f(Nvar, Nvar);
    sigval=malloc(Nvar*sizeof(float));
    iran=randomgenerator();
    InitRandom((RANDOMTYPE)iran);
  }else{
    // clusters k!=clus are copied from previous ones,
    // cluster  clus is randomly split
    if(Round==0){
      if(sort)free(sort);
      sort=Rank_clus(detsig, nc1);
    }
    int clus=sort[Round];
    printf("%d clusters Round %d splitclus= %d\n", ncl, Round, clus);
    int *num=malloc(ncl*sizeof(int));
    for(k=0; k<ncl; k++)num[k]=0;
    for(i=0; i<N; i++){
      for(k=0; k<ncl; k++)w1_ik[i][k]=0;
      int k0=cluster[i];
      if((k0==clus)&&(RandomFloating()<0.5))k0=nc1;
      w1_ik[i][k0]=1; num[k0]++;
    }
    printf("Initial cluster size: ");
    for(k=0; k<=nc1; k++)printf(" %d", num[k]); printf("\n"); 
    free(num);
    // transition probabilities
    if(trans)Compute_trans(trans, w2_ikk1, cluster, ncl, clus, N);
  }


  // Compute initial parameters from weights w1_ik
  Gaussian_parameters_1(tau, mu, sig, w1_ik, x_SamVar, N, Nvar, ncl);

  // Check if the covariance is singular and regularize if not
  int m;
  for(m=0; m<ncl; m++){
    float **sigm=sig[m];
    float det=Determinant(sigm, Nvar);
    if(det<=0.0000001){
      printf("WARNING, ncl=%d ", ncl);
      printf("regularizing singular initial covariance matrix %d tau=%.3g",
	     m, tau[m]);
      f_Diagonalize(Nvar, sigm, sigval, sigvec);
      int n=0;
      for(k=0; k<Nvar; k++){
	if(sigval[k]<lambda0){
	  printf(" Singular component %d, sigma^2=%.5f\n", k, sigval[k]);
	  sigval[k]=lambda0; n++;
	}
      }
      printf("Det= %.6f  %d singular components < %.5f\n", det, n, lambda0);
      for(j=0; j<Nvar; j++){
	for(j1=0; j1<=j; j1++){
	  double sum=0;
	  for(k=0; k<Nvar; k++)sum+=sigval[k]*sigvec[k][j1]*sigvec[k][j];
	  sigm[j][j1]=sum; sigm[j1][j]=sum; 
	}
      }
    }
  }

  return(0);
}

int *Rank_clus(float *detsig, int nc1)
{
  int *sort=malloc(nc1*sizeof(int)), j, k;
  int *sortit=malloc(nc1*sizeof(int));
  for(j=0; j<nc1; j++)sortit[j]=1;
  for(j=0; j<nc1; j++){
    float detmax=-1; int kmax=-1;
    for(k=0; k<nc1; k++){
      if((sortit[k])&&(detsig[k] > detmax)){
	kmax=k; detmax=detsig[k];
      }
    }
    if(kmax < 0){
      printf("ERROR, could not rank previous clusters at rank %d\n", j);
      printf("Cluster to_be_sorted det(sig)\n");
      for(k=0; k<nc1; k++){
	printf("%d %d %f\n", k, sortit[k], detsig[k]);
      }
      exit(8);
    }
    sort[j]=kmax;
    sortit[kmax]=0;
  }
  free(sortit);
  return(sort);
}

void Compute_trans(float **corr, float ***w2_ikk1,
		   int *cluster, int ncl, int clus, int N)
{
  printf("Computing initial transition probabilities\n");
  int nc1=ncl-1, k, k1, i;
  for(k=0; k<ncl; k++){
    for(k1=0; k1<ncl; k1++)corr[k][k1]=0;
  }
  int k0=cluster[0];
  for(i=0; i<N; i++){
    float **w=w2_ikk1[i];
    for(k=0; k<ncl; k++){
      for(k1=0; k1<ncl; k1++)w[k][k1]=0;
    }
    k1=cluster[i];
    if(i==0){w[k1][k1]=1; continue;}
    if((k1==clus)&&(RandomFloating()<0.5))k1=nc1;
    w[k0][k1]=1;
    corr[k0][k1]++;
    k0=k1;
  }
  for(k=0; k<ncl; k++){
    for(k1=0; k1<k; k1++){
      corr[k][k1]=(corr[k][k1]+corr[k1][k])/2;
      corr[k1][k]=corr[k][k1];
    }
  }
  for(k=0; k<ncl; k++){
    long sum=0;
    for(k1=0; k1<ncl; k1++)sum+=corr[k][k1];
    for(k1=0; k1<ncl; k1++)corr[k][k1]/=sum;
  }
}

unsigned long randomgenerator(void){
     
     unsigned long tm;
     time_t seconds;
     
     time(&seconds);
     srand((unsigned)(seconds % 65536));
     do   /* waiting time equal 1 second */
       tm= clock();
     while (tm/CLOCKS_PER_SEC < (unsigned)(1));
     
     return((unsigned long) (rand()));

}

void Cluster_size(int *cluster, int N, int ncl){
  int *num=malloc(ncl*sizeof(int)), i, k;
  for(k=0; k<ncl; k++)num[k]=0;
  for(i=0; i<N; i++)num[cluster[i]]++;
  printf("Cluster sizes: ");
  for(k=0; k<ncl; k++)printf(" %d", num[k]); printf("\n"); 
  free(num);
}
