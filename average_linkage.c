#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "average_linkage.h"
#include "allocate.h"
#define NMIN 20

float Average_linkage_coord(float **Clust_sim, int *cluster, int *Nnew,
			    int N, float **x, int *Nele, int Nvar, float Thr)
{
  int i, k, j, imax, jmax;
  for(i=0; i<N; i++)cluster[i]=i;

  int   *Max_i=malloc(N*sizeof(int));
  float *Max_s=malloc(N*sizeof(float)), s;

  float Max_sim=-1;
  for(i=0; i<N; i++)Max_s[i]=0; 
  for(i=0; i<N; i++){
    Clust_sim[i][i]=1.00;
    for(j=i+1; j<N; j++){
      s=Cosine(x[i], x[j], Nvar);
      Clust_sim[i][j]=s;
      if(s>Max_s[i]){Max_s[i]=s; Max_i[i]=j;}
      if(s>Max_s[j]){Max_s[j]=s; Max_i[j]=i;}
      if(s>Max_sim){Max_sim=s; imax=i; jmax=j;}
    }
  }
  *Nnew=N;
  if(Max_sim < Thr )goto end;
  printf("Average  linkage, %d clusters Thr=%.2f Max_sim=%.2f\n",
	 N, Thr, Max_sim);
  int i1_max, j1_max;
  while(Max_sim>=Thr){
    // Join clusters
    printf("Joining clusters %d %d sim= %.2f N= %d %d\n",
	   imax, jmax, Max_sim, Nele[imax], Nele[jmax]);
    float w1=Nele[imax]/(float)(Nele[imax]+Nele[jmax]), w2=1-w1;
    float *x1=x[imax], *x2=x[jmax];
    for(j=0; j<N; j++)x1[j]=w1*x1[j]+w2*x2[j];
    Nele[imax]+=Nele[jmax]; Nele[jmax]=0;
    int ck=cluster[jmax];
    for(j=0; j<N; j++){
      if(cluster[j]==ck)cluster[j]=imax;
    }
    (*Nnew)--; if((*Nnew)==1)break;

    // Compute similarities
    Max_s[imax]=-1.1; Max_sim=-1;
    for(k=0; k<N; k++){
      if((Nele[k]==0)||(k==imax)||(k==jmax))continue;
      s=Cosine(x[imax], x[k], Nvar);
      if(k<imax){Clust_sim[k][imax]=s;}
      else{Clust_sim[imax][k]=s;}
      if(s>Max_s[imax]){Max_s[imax]=s; Max_i[imax]=k;}

      // Find maximum for k
      if((Max_i[k]==imax)||(Max_i[k]==jmax)){
	Max_s[k]=s; Max_i[k]=imax;
	for(j=0; j<k; j++){
	  if((Nele[j]==0)||(j==imax)||(j==jmax))continue;
	  if(Clust_sim[j][k] > Max_s[k]){
	    Max_s[k]=Clust_sim[j][k]; Max_i[k]=j;
	  }
	}
	for(j=k+1; j<N; j++){
	  if((Nele[j]==0)||(j==imax)||(j==jmax))continue;
	  if(Clust_sim[k][j] > Max_s[k]){
	    Max_s[k]=Clust_sim[k][j]; Max_i[k]=j;
	  }
	}
      }else{
	if(s>Max_s[k]){Max_s[k]=s; Max_i[k]=imax;}
      }
      if(Max_s[k]>Max_sim){
	Max_sim=Max_s[k]; i1_max=k; j1_max=Max_i[k];
      }
    } // End loop of k
    if(i1_max < j1_max){imax=i1_max; jmax=j1_max;}
    else{imax=j1_max; jmax=i1_max;}
  }

  // Reorder cluster indexes
  int cnew=0, cold=0, j1;
  int *newclus=malloc(N*sizeof(int));
  int *oldclus=malloc(N*sizeof(int));
  newclus[0]=0; oldclus[0]=0;
  for(k=0; k<N; k++){
    int ck=cluster[k];
    if(ck > cold){
      cold=ck; cnew++; newclus[ck]=cnew; oldclus[cnew]=ck;
    }
  }
  cnew++;
  if(cnew !=*Nnew){
    printf("ERROR, wrong number of clusters %d %d\n", cnew, *Nnew);
    for(k=0; k<N; k++)printf("c%d: %d  ", k, cluster[k]);
    printf("\n");
    exit(8);
  }

  for(k=0; k<N; k++){
    j=cluster[k]; i=newclus[j];
    if(i!=j){
      cluster[k]=i;
      for(j1=0; j1<Nvar; j1++)x[i][j1]=x[j][j1];
      for(j1=0; j1<i; j1++)Clust_sim[j1][i]=Clust_sim[j1][j];
      for(j1=i+1; j1<N; j1++)Clust_sim[i][j1]=Clust_sim[j][j1];
    }
  }
  for(i=0; i<cnew; i++){
    for(j=i+1; j<N; j++)Clust_sim[j][i]=Clust_sim[i][j];
  }
  free(oldclus); free(newclus);


 end:
  free(Max_i); free(Max_s);
  return(Max_sim);
}    


float Cosine(float *v1, float *v2, int N)
{
  double sum=0, sum1=0, sum2=0; int i;
  float *x=v1, *y=v2;
  for(i=0; i<N; i++){
    sum +=(*x)*(*y);
    sum1+=(*x)*(*x);
    sum2+=(*y)*(*y);
    x++; y++;
  }
  if(sum1 && sum2)sum/=sqrt(sum1*sum2);
  return(sum);
}

int Silouhette(int *cluster, int *numclus, int nclust,
	       float **x_VarSam, int Nsam, int Nvar)
{
  double S=0;
  int k, k1, i, i1, j, j1;
  float **x_SamVar=malloc(Nsam*sizeof(float *));
  for(i=0; i<Nsam; i++){
    x_SamVar[i]=malloc(Nvar*sizeof(float));
    for(j=0; j<Nvar; j++)x_SamVar[i][j]=x_VarSam[j][i];
  }
  double *d_clus=malloc(nclust*sizeof(double));
  int *nswitch=malloc(Nsam*sizeof(int));
  for(i=0; i<Nsam; i++)nswitch[i]=0;
  for(k=0; k<nclust; k++){
    for(i=0; i<Nsam; i++){
      if(cluster[i]!=k)continue;
      for(k1=0; k1<nclust; k1++)d_clus[k1]=0;
      for(i1=0; i1<Nsam; i1++){
	if(i1==i)continue;
	float d, d2=0, *x=x_SamVar[i], *x1=x_SamVar[i1];
	for(j=0; j<Nvar; j++){
	  d=*x-*x1; d2+=d*d; x++; x1++;
	}
	d_clus[cluster[i1]]+=d2;
      } // End loop on other elements
      float Min=100000; int kmin=-1;
      for(k1=0; k1<nclust; k1++){
	if(k1==k){
	  d_clus[k1]/=(numclus[k1]-1);
	}else{
	  d_clus[k1]/=numclus[k1];
	  if(d_clus[k1]<Min){Min=d_clus[k1]; kmin=k1;}
	}
      }
      if(Min<d_clus[k]){
	numclus[k]--; numclus[kmin]++;
	cluster[i]=kmin; nswitch[i]++;	
      }
    }
  }
  int numswitch=0;
  for(i=0; i<Nsam; i++)numswitch+=nswitch[i];
  printf("%d elements switched of cluster in Silouhette\n", numswitch);

  for(k=0; k<nclust; k++){
    if(numclus[k]<NMIN)continue;
  }

  Empty_matrix_f(x_SamVar, Nsam);
  free(nswitch);
  free(d_clus);
}
 
