#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "cluster_score.h"
#include "allocate.h"
#include "expmax.h"

float Cluster_score(int *cluster, int ncl, float **x_VarSam,
		    int Nsam, int Nvar)
{
  int i, j, j2, k;
  int *num=malloc(ncl*sizeof(int));
  int *nfrag=malloc(ncl*sizeof(int));
  double **ave=Allocate_mat2_d(ncl, Nvar);
  double **sig=Allocate_mat2_d(Nvar, Nvar);
  float **sigf=Allocate_mat2_f(Nvar, Nvar);

  double score=0; int new=1;
  // Intracluster score
  for(k=0; k<ncl; k++){
    num[k]=0; nfrag[k]=0;
    for(j=0; j<Nvar; j++){
      for(j2=0; j2<=j; j2++)sig[j][j2]=0;
    }
    double *a=ave[k], *s;
    for(i=0; i<Nsam; i++){
      if((cluster[i]<0)||(cluster[i]>=ncl)){
	printf("ERROR, sample %d has wrong cluster index %d\n", i, cluster[i]);
	continue;
      }
      if(cluster[i]!=k){new=1; continue;}
      num[k]++;
      if(new){nfrag[k]++; new=0;}
      for(j=0; j<Nvar; j++){
	float x=x_VarSam[j][i];
	a[j]+=x; s=sig[j];
	for(j2=0; j2<=j; j2++)s[j2]+=x*x_VarSam[j2][i];
      }
    }
    if(num[k]<1)continue;
    for(j=0; j<Nvar; j++){
      a[j]/=num[k]; s=sig[j];
      for(j2=0; j2<=j; j2++){
	s[j2]=s[j2]/num[k]-a[j]*a[j2];
      }
    }
    for(j=0; j<Nvar; j++){
      for(j2=0; j2<=j; j2++)sigf[j][j2]=sig[j][j2];
      for(j2=j+1; j2<Nvar; j2++)sigf[j2][j]=sigf[j][j2];
    }
    float det=Determinant(sigf, Nvar);
    //score+=(nfrag[k]-1)*log(det);
    if(det>0)score+=(num[k]-1)*log(det);
  }

  // Intercluster score
  if(ncl>1){
    for(j=0; j<Nvar; j++){
      for(j2=0; j2<=j; j2++)sig[j][j2]=0;
    }
    double *avetot=malloc(Nvar*sizeof(double));
    for(j=0; j<Nvar; j++)avetot[j]=0;
    double norm=0, *s;

    for(k=0; k<ncl; k++){
      if(num[k]<1)continue;
      float w=num[k]; norm+=w;
      double *a=ave[k];
      for(j=0; j<Nvar; j++){
	avetot[j]+=w*a[j];
	for(j2=0; j2<=j; j2++)sig[j][j2]+=w*a[j]*a[j2];
      }
    }
    for(j=0; j<Nvar; j++){
      avetot[j]/=norm;
      for(j2=0; j2<=j; j2++){
	sig[j][j2]=sig[j][j2]/norm-avetot[j]*avetot[j2];
      }
    }
    for(j=0; j<Nvar; j++){
      for(j2=0; j2<=j; j2++)sigf[j][j2]=sig[j][j2];
      for(j2=j+1; j2<Nvar; j2++)sigf[j2][j]=sigf[j][j2];
    }
    free(avetot);
    float det=Determinant(sigf, Nvar);
    if(det>0)score+=(ncl-1)*log(det);
  }


  score/=(Nsam-1);

  Empty_matrix_d(ave, ncl);
  Empty_matrix_d(sig, Nvar);
  Empty_matrix_f(sigf, Nvar);
  free(num);
  free(nfrag);
  return(score);
}
