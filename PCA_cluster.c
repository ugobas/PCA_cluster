#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "allocate.h"
#include "expmax.h"
#include "average_linkage.h"
#include "cluster_score.h"

int NORMALIZE_VAR=1; // Normalize variables to compute correlations?
int ORDER=1;         // Use positional order? (HMM)
int nclust_max=12;   // Maximum number of clusters tested
int nclust_min=1;    // Minimum number of clusters tested
float Sim_thr=0.90;  // Threshold for joining clusters
int DIFF_THR=0;      // Minimum number of different variables in cluster
float SIGMA=1;       // Error estimated as SIGMA times standard error of mean
float CHI_THR=1.0;   // Minimum value of Chi to accept a splitting
int NMIN=10;         // Minimum number of elements in a cluster
int NPCA=4;          // Default number of Principal Components for clustering
int COMPACT=0;       // Not used anymore
float SIM_THR;
// Parameters of expectation maximization:
int RMAX=12;          // Maximum number of initial conditions
int IT_MAX=200;      // Maximum number of iterations
float EPS=0.000002;  // Threshold for convergence
float lambda0=0.0001; // Correction to singular matrix 

int Remove_dir(char *nameout, char *filename);
void Get_lcorr(float *lcorr, float **data, int N, int k0, int Nvar);
float Corr_coeff(float *slope, float *offset, float *xx, float *yy, int n);
extern float **Transpose_matrix(float **m, int n1, int n2);

/* Needs name of input file
   Format of the input file: ### with the names of the n variables
   L lines with n variables per line. */

/**************** Input ***************************/
void help(char **argv);
void Get_args(char *file_name, char **file_split, char *name_out,
	      int *N_SHOW, int *Print, int *MATRIX, int *DIFF_THR,
	      float *Sim_thr, int *ORDER, int *nclust_min, int *nclust_max,
	      char *cstop, int argc, char **argv);
float **Read_data(char *file_name, int *N, int *N_data, char ***, char ***,
		  int **select_var);
int Read_var(char *string, int N, int *variables, char **var_names);
int Select_data(float **data, int N, int N_data, char **var_names, int *var);
static int Remove_file_extension(char *nameout, char *ext, char *filename);

/****************** Output *******************/
int Print_PCA(float *, float **, float **, double **,
	      char **, char **, int, int, char *, double *, double *,
	      int *N_SHOW, int Print);
int Print_results(int **clust_store, float **P_store, int *numclus,
		  float **ave, float **err, float *lcorr,
		  float **data, float **x_VarSam, char **names,
		  int Nsam, int Nvar, int NPCA, float **Clust_sim,
		  float Sim_thr, float **Clust_sim_bin, float SIGMA,
		  FILE *filelog, char *txt, char *nameout, int Print,
		  char *file_split, int ncl, int ncl_min, int ncl_max);
static int Print_parameters(int *numcl,
			    float *tau, float **corr,
			    float **mu, float ***sig,
			    int ncl, int Nvar, FILE *file);
static int Print_propensity(float *tau, float **corr, int ncl, FILE *filelog);

char **Set_names(char **var_names, int Nvar, int NPCA);
void Print_selected_data(float **data, int Nvar, int Nsam,
			 char **var_names, char *file_name);
void Print_clusters(int **clust_store, float **P_store, int Nsam,
		    char *file_name, int ncl_min, int ncl_max, int ncl);
int Print_profile(int *numclus, float **ave, float **dev, char **var_names,
		  int nclust, int Nsam, int Nvar, char *nameout, char *open);
int Print_sim(float **Clust_sim, int nclust, int Nvar, char *nameout,
	      char *open, char *what);
int Split_file(char *file_split, int nclust, int *cluster, int Nsam);

/********************** PCA *************************/
extern void Diagonalize(int N, double **MATRIX,
			float *eigen_values, float **eigen_vector);
extern void f_Diagonalize(int N, float **MATRIX,
			float *eigen_values, float **eigen_vector);
int Correlation(double **C_matrix, float **data, int Nvar, int Nsam);
int Compute_corr(double **C_matrix, double *average, double *deviation,
		 float **data_norm, float **data, int Nvar, int Nsam);
int Compute_PC(float **x_VarSam, float **PCA_load, float *eigen_values,
	       float **data, int num, int Nvar, int NPC);

/********************** Computations *************************/
float AIC(float lik, int Npar, int Nsam);
float BIC(float lik, int Npar, int Nsam);
void Cluster_statistics(int *numclus, float **ave, float **dev, int nclust,
			int *cluster, float **data, int Nsam, int j0, int Nvar,
			float *tau);
float Cosine(float *v1, float *v2, int N);
float Chi2(float *ave1, float *ave2, float *dev1, float *dev2, int N);
float Similarity(float *ave1, float *ave2, float *dev1, float *dev2,
		 int N, float SIGMA);

/********************** Clustering *************************/
int Select_points(float **data_select, float **data, int *index,
		  int *cluster, int clust, int Nsam, int Nvar);


/******************* Main routine *******************/

int main(int argc, char **argv)
{
  int Print=0, MATRIX=0;

  float **data;
  char **var_names, **ind_names, file_name[200], name[200], cstop='b';
  int Nvar=0, Nsam, i, j, *select_var;
  char nameout[300], nameout2[300], *file_split=NULL;

  Get_args(file_name, &file_split, nameout, &NPCA, &Print, &MATRIX, &DIFF_THR,
	   &Sim_thr, &ORDER, &nclust_min, &nclust_max, &cstop, argc, argv);
  data=Read_data(file_name, &Nvar, &Nsam, &var_names, &ind_names, &select_var);
  Nvar=Select_data(data, Nvar, Nsam, var_names, select_var);
  if(MATRIX)Print_selected_data(data, Nvar, Nsam, var_names, file_name);
  SIM_THR=(Nvar-DIFF_THR)/((float)Nvar)+0.0001;

  char ext[100];
  Remove_dir(nameout, file_name);
  Remove_file_extension(nameout, ext, nameout);

  // Allocate arrays
  double *average=malloc(Nvar*sizeof(double));
  double *deviation=malloc(Nvar*sizeof(double));
  double **corr_norm=Allocate_mat2_d(Nvar, Nvar);
  float **data_norm=Allocate_mat2_f(Nvar, Nsam);
  Compute_corr(corr_norm, average, deviation, data_norm, data, Nvar, Nsam);
  printf("%d elements, %d selected variables:\n", Nsam, Nvar);
  float norm=(2.0)/sqrt(Nsam);
  for(i=0; i<Nvar; i++){
    printf("%9.3g\t%.1g\t%s\n", average[i], norm*deviation[i], var_names[i]);
  }

  if(NPCA > Nvar)NPCA=Nvar;
  double **correlation=Allocate_mat2_d(Nvar, Nvar), logdet=0;
  float *eigen_values=malloc(Nvar*sizeof(float));
  float **PCA_load=Allocate_mat2_f(Nvar, Nvar);
  float **x_VarSam=Allocate_mat2_f(NPCA, Nsam);
  // Correlation of unnormalized variables
  Correlation(correlation, data, Nvar, Nsam);
  double **corr_dat; float **dat;
  if(NORMALIZE_VAR){
    // Diagonalize correlation of normalized variables
    corr_dat=corr_norm; dat=data_norm;
  }else{
    // Diagonalize correlation of unnormalized variables
    corr_dat=correlation; dat=data;
  }
  Diagonalize(Nvar, corr_dat, eigen_values, PCA_load);
  Compute_PC(x_VarSam, PCA_load, eigen_values, dat, Nsam, Nvar, NPCA);
  Print_PCA(eigen_values, PCA_load, dat, corr_norm,
	    var_names, ind_names, Nvar, Nsam, nameout,
	    average, deviation, &NPCA, Print);
  
  // Correlation length of data
  float *lcorr=malloc((Nvar+NPCA)*sizeof(float));
  Get_lcorr(lcorr, data, Nsam, 0, Nvar);
  Get_lcorr(lcorr, x_VarSam, Nsam, Nvar, NPCA);
  float lcmax=1, Nsam_eff;
  for(i=0; i<NPCA; i++)if(lcorr[Nvar+i]>lcmax)lcmax=lcorr[Nvar+i];
  printf("Maximum correlation length of PCs: %.0f\n", lcmax);
  if((lcmax > 0)&&(lcmax < Nsam)){Nsam_eff=Nsam/lcmax;}
  else{Nsam_eff=Nsam/5;}
  printf("Number of effective independent samples= %.2f\n", Nsam_eff);


  // Output
  for(i=0; i<Nvar; i++)logdet+=log(eigen_values[i]);
  printf("-log(det(correlation))= %.4f\n", -logdet);
  sprintf(nameout, "%s_n%d", nameout, NPCA);
  if(ORDER)sprintf(nameout, "%s_HMM", nameout);
  sprintf(nameout2, "%s.log", nameout);
  FILE *filelog=fopen(nameout2, "w");
  sprintf(nameout2, "%s_bic.dat", nameout);
  FILE *filebic=fopen(nameout2, "w");
  fprintf(filebic, "# Correlation length of PCs: ");
  for(i=0; i<NPCA; i++)fprintf(filebic, "  %d = %.0f", i+1, lcorr[Nvar+i]);
  fprintf(filebic, "\n");
  fprintf(filebic, "# Maximum correlation length of PCs: %.0f\n", lcmax);
  fprintf(filebic, "# Number of effective independent samples= %.0f\n",
	  Nsam_eff);
  fprintf(filebic, "# log(det(correlation))= %.4f\n", logdet);
  fprintf(filebic, "# PCA_cluster %s %d PC %d samples\n", nameout, NPCA, Nsam);
  fprintf(filebic,
	  "# nclus Dlikelihood DAIC DBIC Cluster_score Mean_square_cluster_separation max_Cosine\n");
  char txt[10000];
  sprintf(txt, "Clustering %d PCs\n", NPCA);
  sprintf(txt, "%sVariance explained by the PCs that are used:\n", txt);
  for(i=0; i<NPCA; i++)sprintf(txt,"%s%2d %.2f\n", txt,i,eigen_values[i]);
  sprintf(txt, "%sVariance explained by the PCs that are not used:\n", txt);
  for(i=NPCA; i<Nvar; i++)sprintf(txt, "%s%2d %.2f\n", txt,i,eigen_values[i]);
  printf("%s", txt); fprintf(filelog, "%s", txt);

  // Printing
  char **names=Set_names(var_names, Nvar, NPCA);
  float **Clust_sim=Allocate_mat2_f(nclust_max, nclust_max);
  float **Clust_sim_bin=Allocate_mat2_f(nclust_max, nclust_max);

  // Clustering
  int cluster[Nsam]; for(i=0; i<Nsam; i++)cluster[i]=0;
  int ncl=1, ncl_opt=1, clust, k;
  float lik, aic, bic, cscore, AIC_MIN, BIC_MIN, C_MIN, CHI_MAX;
  float lik1, aic1, bic1, g1;
  int numclus[nclust_max+1];
  float **ave=Allocate_mat2_f((nclust_max+1), Nvar+NPCA);
  float **err=Allocate_mat2_f((nclust_max+1), Nvar+NPCA);

  float *tau=malloc(nclust_max*sizeof(float));
  float **corr=Allocate_mat2_f(nclust_max, nclust_max);
  float **mu=Allocate_mat2_f(nclust_max, NPCA);
  float ***sig=malloc(nclust_max*sizeof(float **));
  for(k=0; k<nclust_max; k++)sig[k]=Allocate_mat2_f(NPCA, NPCA);
  float Pstate[Nsam];
  int join[nclust_max];

  // Storing_results
  int nclust_store=nclust_max+1;
  int *clust_store[nclust_store];
  float *Pstate_store[nclust_store];
  for(i=1; i<nclust_store; i++){
    clust_store[i]=malloc(Nsam*sizeof(int));
    Pstate_store[i]=malloc(Nsam*sizeof(float));
  }

  float norm_lik1=Nsam_eff/Nsam, norm_lik2=Nsam_eff*NPCA;
  for(ncl=nclust_min; ncl<=nclust_max; ncl++){
    int split=1, first; char open[4];
    if(ncl==nclust_min){
      strcpy(open,"w"); first=1;
    }else{
      strcpy(open,"a"); first=0;
    }
    if(ORDER){
      lik=Expectation_maximization_order(cluster, Pstate,
					 ncl, tau, corr, mu, sig,
					 x_VarSam, Nsam, NPCA, first);
    }else{
      lik=Expectation_maximization(cluster, Pstate,
				   ncl, tau, mu, sig,
				   x_VarSam, Nsam, NPCA, first);
    }

    // Compute likelihood and other scores
    lik*=norm_lik1;
    int N_para = ncl*(NPCA + NPCA*(NPCA+1)/2 +1)-1; // mu, sigma, tau
    if(ORDER)N_para+=(ncl*(ncl-1));
    aic=AIC(lik, N_para, Nsam_eff);
    bic=BIC(lik, N_para, Nsam_eff);
    lik/=norm_lik2; aic/=norm_lik2; bic/=norm_lik2;
    float cscore=Cluster_score(cluster, ncl, x_VarSam, Nsam, NPCA);
    if(ncl==nclust_min){lik1=lik; aic1=aic; bic1=bic; g1=cscore;}
    sprintf(txt,
	    "%d clusters: lik= %.4f AIC=%.4f BIC= %.4f (/NPC) score=%.3f\n",
	    ncl, lik, aic, bic, cscore);
    printf("%s", txt); fprintf(filelog, "%s", txt);
    int N_clus=ncl;

    // Separation score
    Cluster_statistics(numclus,ave,err,ncl,cluster,dat,Nsam,0,Nvar,lcorr);
    double chi, chisum=0, simsum=0, sum=0, q, sim_max=-1;
    //if(ncl==1)goto print_bic;
    float w[ncl], ww;
    float Chi_Thr=CHI_THR/sqrt(((float)Nsam/(float)ncl));
    for(i=0; i<ncl; i++)join[i]=-1;
    for(i=0; i<ncl; i++){
      //if(numclus[i]<NMIN)continue;
      int sim=0;
      w[i]=(float)numclus[i]/Nsam;
      for(j=i+1; j<ncl; j++){
	//if(numclus[j]<NMIN)continue;
	chi=Chi2(ave[i], ave[j], err[i], err[j], Nvar)/Nvar;
	ww=w[i]*w[j]; chisum+=ww*chi; sum+=ww;
	if(chi < Chi_Thr){
	  sprintf(txt,"Too similar clusters %d %d, chi= %.2f\n", i,j,chi);
	  printf("%s", txt); fprintf(filelog, "%s", txt);
	  split=0; sim=1;
	}
	q=Similarity(ave[i], ave[j], err[i], err[j], Nvar, SIGMA);
	if(q > SIM_THR){
	  sprintf(txt, "Too similar clusters %d %d, differences= %.1f ",
		 i, j, (1-q)*Nvar);
	  sprintf(txt, "%sThreshold= %d\n", txt, DIFF_THR);
	  printf("%s", txt); fprintf(filelog, "%s", txt);
	  split=0; sim=1;
	}
	q=Cosine(ave[i], ave[j], Nvar);
	simsum+=ww*q; if(q > sim_max)sim_max=q;
      }
      if(sim)N_clus--;
    }
    if(sum){chisum/=sum; simsum/=sum;}
    sprintf(txt,
	    "Clusters after joining: %d Discriminative power, chi2= %.3f\n\n",
	    N_clus,chisum);
    printf("%s", txt); fprintf(filelog, "%s", txt);

    // Optimal number of clusters
    // Optimal scores
    int s_bic=0, s_aic=0, s_score=0, s_sep=0;
    if(ncl==nclust_min){
      AIC_MIN=aic; BIC_MIN=bic; C_MIN=cscore; CHI_MAX=chisum;
      s_bic=1, s_aic=1, s_score=1, s_sep=1;
    }else{
      if(aic < AIC_MIN){AIC_MIN=aic; s_aic=1;}
      if(bic < BIC_MIN){BIC_MIN=bic; s_bic=1;}
      if(cscore <C_MIN){C_MIN=cscore; s_score=1;}
      if(chisum > CHI_MAX){CHI_MAX=chisum; s_sep=1;}
    }

    if(((cstop=='a')&&(s_aic))||
       ((cstop=='b')&&(s_bic))||
       ((cstop=='c')&&(s_score))||
       ((cstop=='s')&&(s_sep))){
      ncl_opt=ncl;
    }else { // Stop splitting clusters
      //printf("Stopping clustering, optimal number= %d\n", ncl_opt);
      split=0;
    }

    // Test if small clusters
    int n2_clus=N_clus;
    for(i=0; i<ncl; i++)if(numclus[i]<NMIN){n2_clus--; split=0;}
    if(n2_clus < N_clus){
      sprintf(txt, "WARNING, %d clusters have fewer than %d elements\n",
	      N_clus-n2_clus, NMIN);
      printf("%s", txt); fprintf(filelog, "%s", txt);
    }

  print_bic:
    fprintf(filebic, "%2d %7.4f %7.4f %7.4f %7.4f %.3f %6.3f\n", N_clus,
	    lik-lik1, aic-aic1, bic-bic1, cscore-g1, chisum, sim_max);
    fflush(filebic);
    /*if(split==0){
      sprintf(txt, "Stopping clustering, %d clusters\n", ncl-1);
      printf("%s", txt); fprintf(filelog, "%s", txt);
      break;
    }
    */
    if(ORDER){
      Print_parameters(numclus, tau, corr, mu, sig, ncl, NPCA, filelog);
      Print_propensity(tau, corr, ncl, filelog);
    }
    int *cl=clust_store[ncl]; float *P=Pstate_store[ncl];
    for(i=0; i<Nsam; i++){cl[i]=cluster[i]; P[i]=Pstate[i];}
    Print_results(clust_store, Pstate_store, numclus, ave, err, lcorr,
		  dat, x_VarSam, names, Nsam, Nvar, NPCA, Clust_sim, Sim_thr,
		  Clust_sim_bin, SIGMA, filelog, txt, nameout, Print,
		  file_split, ncl, nclust_min, ncl);

  }
  fclose(filebic);

  printf("Optimal number of clusters: %d\n", ncl_opt);
  if(ncl_opt!=nclust_max){
    Print_results(clust_store, Pstate_store, numclus, ave, err, lcorr,
		  dat, x_VarSam, names, Nsam, Nvar, NPCA, Clust_sim, Sim_thr,
		  Clust_sim_bin, SIGMA, filelog, txt, nameout, Print,
		  file_split, ncl_opt, nclust_min, nclust_max);
  }

  printf("Parameters used: %d variables, %d principal components, %d samples\n",
	 Nvar, NPCA, Nsam);
  printf("Minimum number of different variables in clusters: %d\n", DIFF_THR);
  printf("Minimum number of elements in a cluster: %d\n", NMIN);
  printf("Error estimated as %.1f times the standard deviation\n", SIGMA);
  printf("Threshold similarity: %.3f\n", Sim_thr);

  Empty_Gaussian_parameters(tau,mu,sig, NPCA, nclust_max);
  Empty_matrix_f(corr, nclust_max);
  return(0);
}

int Print_results(int **clust_store, float **P_store, int *numclus,
		  float **ave, float **err, float *lcorr,
		  float **data, float **x_VarSam, char **names,
		  int Nsam, int Nvar, int NPCA, float **Clust_sim,
		  float Sim_thr, float **Clust_sim_bin, float SIGMA,
		  FILE *filelog, char *txt, char *nameout, int Print,
		  char *file_split, int ncl_act, int ncl_min, int ncl_max)
{
  if(ncl_act==1)return(0);
  int ncl=ncl_act, i, j, *cluster=clust_store[ncl];
  printf("Printing results for %d clusters\n", ncl);

  //if(ncl != ncl_max){ // Optimal number of clusters, join similar clusters
  // Join similar clusters

  // Compute average observables per cluster (original and PCs)
  Cluster_statistics(numclus, ave, err, ncl, cluster,
		     data, Nsam, 0, Nvar, lcorr);
  Cluster_statistics(numclus, ave, err, ncl, cluster,
		     x_VarSam, Nsam, Nvar, NPCA, lcorr);

  // Join similar clusters through average linkage
  int clus_clus[ncl], nnew=ncl;
  float Max_sim=
    Average_linkage_coord(Clust_sim, clus_clus, &nnew,
			  ncl, ave, numclus, Nvar, Sim_thr);
  if(nnew != ncl){
    sprintf(txt,"Average linkage, joining %d clusters into %d\n\n",ncl,nnew);
    printf("%s", txt); fprintf(filelog, "%s", txt);
    ncl=nnew; for(i=0; i<Nsam; i++)cluster[i]=clus_clus[cluster[i]];
    Cluster_statistics(numclus, ave, err, ncl, cluster,
		       data, Nsam, 0, Nvar, lcorr);
    Cluster_statistics(numclus, ave, err, ncl, cluster,
		       x_VarSam, Nsam, Nvar, NPCA, lcorr);
  }

  // Rename clusters according to the first PC
  int O_all=Nvar+NPCA, k;
  float score[ncl]; for(k=0; k<ncl; k++)score[k]=ave[k][Nvar];
  Sort(clus_clus, score, ncl);
  for(i=0; i<Nsam; i++)cluster[i]=clus_clus[cluster[i]];
  for(i=0; i<O_all; i++){
    float ave_tmp[ncl], err_tmp[ncl];
    for(k=0; k<ncl; k++){
      ave_tmp[clus_clus[k]]=ave[k][i];
      err_tmp[clus_clus[k]]=err[k][i];
    }
    for(k=0; k<ncl; k++){
      ave[k][i]=ave_tmp[clus_clus[k]];
      err[k][i]=err_tmp[clus_clus[k]];
    }
  }

  //}

  // Silouhette
  // Silouhette(cluster, numclus, ncl, data, Nsam, Nvar);

  // Compute similarity matrices
  for(i=0; i<ncl; i++){
    Clust_sim[i][i]=1.0;
    Clust_sim_bin[i][i]=1.0;
    for(j=0; j<i; j++){
      Clust_sim[i][j]=Cosine(ave[i], ave[j], Nvar);
      Clust_sim_bin[i][j]=Similarity(ave[i],ave[j],err[i],err[j],Nvar,SIGMA);
      Clust_sim_bin[j][i]=Clust_sim_bin[i][j];
      Clust_sim[j][i]=Clust_sim[i][j];
    }
  }

  // Printing
  Print_clusters(clust_store, P_store, Nsam, nameout, ncl, ncl_min, ncl_max);
  
  char open[4];
  if(ncl==2 || ncl==ncl_min){strcpy(open,"w");}else{strcpy(open,"a");}
  Print_profile(numclus, ave, err, names, ncl, Nsam, Nvar+NPCA, nameout,open);
  Print_sim(Clust_sim, ncl, Nvar, nameout, open,"Cosine(ave_i,ave_j)");
  // Continuous
  Print_sim(Clust_sim_bin, ncl, Nvar, nameout, "a", "BinaryCos(ave_i,ave_j)");
  // Binary

  if(Print){
    char namesplit[400];
    sprintf(namesplit, "%s_split%d.dat", nameout, ncl);
    Print_cluster_PC(x_VarSam, cluster, ncl, Nsam, NPCA, namesplit);
  }
  if(file_split)Split_file(file_split, ncl, cluster, Nsam);
  printf("\n");
  return(0);
}

float Similarity(float *ave1, float *ave2, float *dev1, float *dev2,
		 int N, float SIGMA)
{
  int j, k1, k2;
  double q=0, n1=0, n2=0, dev;
  for(j=0; j<N; j++){
    dev=dev1[j]; if(dev2[j]>dev)dev=dev2[j];
    dev*=SIGMA;
    if(ave1[j]>dev){k1=1;}
    else if(ave1[j]<-dev){k1=-1;}
    else{k1=0;}
    if(ave2[j]>dev){k2=1;}
    else if(ave2[j]<-dev){k2=-1;}
    else{k2=0;}
    //q+=k1*k2; n1+=k1*k1; n2+=k2*k2;
    if((k1==k2))q++;
  }
  //q/=sqrt(n1*n2);
  return(q/N);
}

float Chi2(float *ave1, float *ave2, float *dev1, float *dev2, int N){
  int j; double chi=0;
  for(j=0; j<N; j++){
    if((dev1[j]>0)&&(dev2[j]>0)){
      float x=ave1[j]-ave2[j];
      chi+=x*x;
      //chi+=(x*x/(dev1[j]*dev2[j]));
    }
  }
  return(chi);
}

void Get_args(char *file_name, char **file_split, char *name_out,
	      int *N_SHOW, int *Print, int *MATRIX, int *DIFF_THR,
	      float *Sim_thr, int *ORDER, int *nclust_min, int *nclust_max,
	      char *cstop, int argc, char **argv)
{
  int i;
  if(argc<2)help(argv);
  for(i=1; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0){
      help(argv);
    }else if(i==1){
      strcpy(file_name, argv[1]);
    }else if(strncmp(argv[i], "-stop", 5)==0){
      char stop[80];
      i++; sscanf(argv[i], "%s", stop);
      if(strncmp(stop, "aic", 3)==0){
	*cstop='a'; printf("Changing stopping criterion to aic\n");
      }else if(strncmp(stop, "bic", 3)==0){
	*cstop='b'; printf("Changing stopping criterion to bic\n");
      }else if(strncmp(stop, "score", 5)==0){
	*cstop='c'; printf("Changing stopping criterion to score\n");
      }else if(strncmp(stop, "sep", 3)==0){
	*cstop='s'; printf("Changing stopping criterion to separation\n");
      }else{
	printf("WARNING, %s is not an allowed stopping criterion\n", stop);
	printf("Using default: bic\n");
      }
    }else if(strncmp(argv[i], "-smin", 5)==0){
      i++; sscanf(argv[i], "%d", nclust_min);
      printf("Minimum number of tested states: %d\n", *nclust_min);
    }else if(strncmp(argv[i], "-smax", 5)==0){
      i++; sscanf(argv[i], "%d", nclust_max);
      printf("Maximum number of tested states: %d\n", *nclust_max);
    }else if(strncmp(argv[i], "-diff", 5)==0){
      i++; sscanf(argv[i], "%d", DIFF_THR);
    }else if(strncmp(argv[i], "-print", 6)==0){
      *Print=1;
    }else if(strncmp(argv[i], "-matrix", 6)==0){
      *MATRIX=1;
    }else if(strncmp(argv[i], "-noorder", 6)==0){
      *ORDER=0;
      printf("Do not take into account positional correlations (HMM)\n");
    }else if(strncmp(argv[i], "-split", 6)==0){
      i++; *file_split=malloc(400*sizeof(char));
      sscanf(argv[i], "%s", *file_split);
    }else if(strncmp(argv[i], "-sim", 4)==0){
      i++; sscanf(argv[i], "%f", Sim_thr);
      printf("Similarity threshold set to %.3f\n", *Sim_thr);
    }else if(strncmp(argv[i], "-n", 2)==0){
      i++; sscanf(argv[i], "%d", N_SHOW);
      if(*N_SHOW >0){
	printf("Number of modes: %d\n", *N_SHOW);
      }else{
	printf("WARNING, %d not allowed as number of modes",*N_SHOW);
	printf("set to default %d\n", NPCA);
	*N_SHOW=NPCA;
      }
    }else{
      printf("WARNING, i=%d unknown option %s\n", i, argv[i]);
    }
  }
  if(file_name[0]==0)help(argv);
  if(*nclust_min > *nclust_max){
    printf("WARNING, minimum number of states %d larger than maximum %d\n",
	   *nclust_min, *nclust_max);
    *nclust_min=1;
    printf("Changing minimum to %d\n", *nclust_min);
  }
}

void help(char **argv){
  printf("\nProgram %s, author Ugo Bastolla.\n", argv[0]);
  printf("Usage: %s <input_file> (mandatory)\n", argv[0]);
  printf("Input: matrix with N rows (samples) and M columns (variables).\n");
  printf("First line must contain variable names:\n");
  printf("  ### name1 name2 name3...\n");
  printf("Second line selects variables (1=selected, 0=omitted):\n");
  printf("  ### 1 0 1 ...\n");
  printf("Options:\n");
  printf("   -h help\n");
  printf("   -smin <Minimum number of states> default %d\n", nclust_min);
  printf("   -smax <Maximum number of states> default %d\n", nclust_max);
  printf("   -stop <stopping criterion> Allowed: aic, bic, sep, score\n");
  printf("   -n <number of PCA used for clustering> (default %d)\n", NPCA);
  printf("   -noorder: Do not take into account positional correlations (HMM)\n");
  printf("   -sim <Maximum similarity of two clusters>, default %.3f\n",
	 Sim_thr);
  printf("   -diff <Minimum number of different variables in two clusters>");
  printf(" (default %d)\n", DIFF_THR);
  printf("   -split <file> file with N lines (as samples) to be splitted\n");
  printf("   -print  Print n principal components of all elements\n");
  printf("   -matrix Print matrix with selected variables\n");
  printf("\n");
  exit(8);
}

float ** Read_data(char *file_name, int *N, int *N_data,
		   char ***var_names, char ***ind_names,
		   int **select_var)
{
  FILE *file_in=fopen(file_name, "r");
  float **data;
  char string[2000]; int i, m, m_data=0, read_names=0, read_var=0;
  *N=0; *N_data=0;  *ind_names=NULL;

  if(file_in==NULL){
    printf("Error, file %s not found\n", file_name); exit(8);
  }

  while(fgets(string, sizeof(string), file_in)!=0){
    if(string[0]=='#')continue;
    if(*N_data==0)*N=Count_columns(string);
    (*N_data)++;
  }
  fclose(file_in);
  printf(" %d columns of %d lines found in %s\n", *N, *N_data, file_name);

  /* Prepare variables */
  data=malloc((*N)*sizeof(float *));
  for(i=0; i< (*N); i++)data[i]=malloc((*N_data)*sizeof(float));
  *var_names=malloc((*N)*sizeof(char *));
  for(i=0; i< (*N); i++){
    (*var_names)[i]=malloc(500*sizeof(char));
    sprintf((*var_names)[i], "var%d\0", i+1); 
  }
  *select_var=malloc((*N)*sizeof(int));
  for(i=0; i<(*N); i++)(*select_var)[i]=1;

  i=0; m_data=0;
  file_in=fopen(file_name, "r");
  while(fgets(string, sizeof(string), file_in)!=0){
    if(string[0]=='#'){
      if(strncmp(string, "###", 3)==0){
	if(read_names==0){
	  Read_names(string, *N, *var_names); read_names=1;
	}else if(read_var==0){
	  Read_var(string, *N, *select_var, *var_names); read_var=1; 
	}
      }else if(m_data>0){
	/*if(*ind_names==NULL){
	  *ind_names=malloc((*N_data)*sizeof(char *));
	  for(i=0; i< (*N_data); i++){
	    (*ind_names)[i]=malloc(100*sizeof(char));
	    sprintf((*ind_names)[i], "I%d\0", i+1);
	  }
	}
	sscanf(string+1, "%s", (*ind_names)[m_data-1]);*/
      }
    }else{
      m=Read_columns(string, *N, data, m_data);
      if(m!=*N){
	printf("Error in line %d, %d columns found\n", i+1, m); m_data--;
      }
      m_data++; i++;
    }
  }
  (*N_data)=m_data;
  return(data);
}


int Count_columns(char *string){
  int m=0, read=0;
  char *ptr=string;

  while(*ptr!='\n'){
    if((read==0)&&(*ptr!=' ')&&(*ptr!='\t')){
      read=1; m++;
    }else if(((*ptr==' ')||(*ptr=='\t'))&&(read==1)){
      read=0;
    }
    ptr++;
  }
  return(m);
}


int Read_columns(char *string, int N, float **data, int m_data){
  int m=0, read=0;
  char *ptr=string;

  while(*ptr!='\n'){
    if((read==0)&&(*ptr!=' ')&&(*ptr!='\t')){
      read=1;
      if(m>=N)return(m+1);
      sscanf(ptr, "%f", &(data[m][m_data]));
      m++;
    }else if(((*ptr==' ')||(*ptr=='\t'))&&(read==1)){
      read=0;
    }
    ptr++;
  }
  return(m);
}

int Read_names(char *string, int N, char **var_names){
  int m=0, read=0;
  char *ptr=string+3;

  while(*ptr!='\n'){
    if((read==0)&&(*ptr!=' ')&&(*ptr!='\t')){
      read=1;
      if(m>=N)return(m+1);
      sscanf(ptr, "%s", var_names[m]);
      m++;
    }else if(((*ptr==' ')||(*ptr=='\t'))&&(read==1)){
      read=0;
    }
    ptr++;
  }
  return(m);
}

int Read_var(char *string, int N, int *variables, char **var_names)
{
  int m=0, read=0;
  char *ptr=string+3;

  while(*ptr!='\n'){
    if((read==0)&&(*ptr!=' ')&&(*ptr!='\t')){
      read=1;
      if(m>=N)return(m+1);
      sscanf(ptr, "%d", &variables[m]);
      m++;
    }else if(((*ptr==' ')||(*ptr=='\t'))&&(read==1)){
      read=0;
    }
    ptr++;
  }
  return(m);

}

int Print_PCA(float *eigen_values, float **eigen_vector, float **data_norm,
	      double **correlation, char **var_names, char **ind_names,
	      int N, int N_data, char *name_out, double *average,
	      double *deviation, int *N_SHOW, int Print)
{
  if(*N_SHOW>=N)*N_SHOW=N-1;

  char filename[400]; sprintf(filename, "%s.pca", name_out);
  FILE *file_out=fopen(filename, "w");
  printf("Writing %s\n", filename);
  fprintf(file_out, "# %d variables, %d elements %d PC\n",
	  N, N_data, *N_SHOW);

  int v1, v2, e1, i;
  float r1, r2, *mean=malloc(N*sizeof(float));
  fprintf(file_out, "# Correlation matrix\n");
  for(v1=0; v1<N; v1++){
    for(v2=v1+1; v2<N; v2++){
      //fprintf(file_out, "# %6.3f ", correlation[v1][v2]);
      //fprintf(file_out, "%12s %12s\n", var_names[v1], var_names[v2]);
      fprintf(file_out, "%s\t%s\t%.3f\n",
	      var_names[v1], var_names[v2],
	      correlation[v1][v2]);
    }
  }
  if(eigen_values==NULL){fclose(file_out); return(0);}


  double sum=0; for(e1=0; e1<N; e1++)sum+=eigen_values[e1];
  double tot=0; int N_sig=0;
  fprintf(file_out, "# Spectrum: n lambda_n cumulative\n");
  for(e1=0; e1<N; e1++){
    tot+=eigen_values[e1];
    double e=0; for(i=0; i<N; i++)e+=eigen_vector[i][e1]; e/=N;
    if(e<0){
      for(i=0; i<N; i++)eigen_vector[i][e1]=-eigen_vector[i][e1]; e=-e;
    }
    mean[e1]=e;
    float s=(eigen_values[e1]/sum)*N;
    if(Print)fprintf(file_out, "# ");
    fprintf(file_out, "%2d\t%6.3f\t%.3f\n",e1, s, tot/sum);
    if(s>1.0)N_sig++;
  }
  if(*N_SHOW==0)*N_SHOW=N_sig;

  fprintf(file_out, "# %d variables, %d elements\n", N, N_data);
  fprintf(file_out, "#      Ave       s.e.       sigma  ");
  for(e1=0; e1<*N_SHOW; e1++)fprintf(file_out, " PC%d   ", e1);
  fprintf(file_out, "sum(PC)^2  var.\n");
  for(v1=0; v1<N; v1++){
    if(Print)fprintf(file_out, "# ");
    fprintf(file_out, "%.3f\t%.3f\t%.3f",
	    average[v1], deviation[v1]/sqrt(N), deviation[v1]);
    tot=0;
    for(e1=0; e1<*N_SHOW; e1++){
      if(e1>=N)break;
      r1= eigen_vector[v1][e1]*sqrt(eigen_values[e1]);
      //r1= eigen_vector[v1][e1]*sqrt(eigen_values[e1])/deviation[v1];
      fprintf(file_out, "\t%.3f", r1);
      tot+=r1*r1;
    }
    fprintf(file_out, "\t%.3f\t%s\n",
	    tot/(deviation[v1]*deviation[v1]), var_names[v1]);
    //fprintf(file_out, "\t%.3f\t%s\n", tot, var_names[v1]);
  }

  if(Print){
    fprintf(file_out, "# Components 1...%d  Individual\n", *N_SHOW);
    for(i=0; i<N_data; i++){
      for(e1=0; e1<*N_SHOW; e1++){     
	if(e1>=N)break; r1=0;
	for(v1=0; v1<N; v1++)
	  r1+=eigen_vector[v1][e1]*(data_norm[v1][i]);
	fprintf(file_out, "%6.3f ", r1);
      }
      if(ind_names)fprintf(file_out, "  %s", ind_names[i]);
      fprintf(file_out, "\n");
    }
  }

  fclose(file_out);
  return(0);
}


int Compute_corr(double **C_matrix, double *average, double *deviation,
		 float **data_norm, float **data, int Nvar, int Nsam)
{
  double corr=0, ave, dev;
  float *dat1, *dat2;
  int k1, k2, i;

  for(k1=0; k1<Nvar; k1++){
    ave=0; dev=0;
    dat1=data[k1]; C_matrix[k1][k1]=1.;
    for(i=0; i<Nsam; i++){
      ave+=dat1[i]; dev+=dat1[i]*dat1[i];
    }   
    ave/=Nsam; dev=(dev/Nsam)-(ave*ave); dev=sqrt(dev);
    for(i=0; i<Nsam; i++){
      data_norm[k1][i]=(dat1[i]-ave)/dev;
    }
    average[k1]=ave; deviation[k1]=dev;

    for(k2=0; k2< k1; k2++){
      corr=0; dat2=data[k2];
      for(i=0; i<Nsam; i++){
	corr+=dat1[i]*dat2[i];
      }
      corr=(corr/Nsam-(ave*average[k2]))/(dev*deviation[k2]);
      C_matrix[k1][k2]=corr; C_matrix[k2][k1]=corr;
    }
  }
  return(0);
}


int Select_data(float **data, int N, int N_data, char **var_names, int *var)
{
  int i, j, k, M=0;
  char **data_tmp, **name_tmp;
  for(i=0; i<N; i++)if(var[i])M++;
  printf("Selecting %d variables out of %d\n", M, N);
  if(M==N)return(M);
  i=0;
  for(k=0; k<M; k++){
    while(var[i]==0){
      printf("Omitting: %s\n", var_names[i]);
      i++;
    }
    if(i!=k){
      for(j=0; j<N_data; j++)data[k][j]=data[i][j];
      strcpy(var_names[k], var_names[i]);
    }
    i++;
  }
  for(j=i; j<N; j++)
    if(var[j]==0)printf("Omitting: %s\n", var_names[j]);
  return(M);
}

void Print_selected_data(float **data, int Nvar, int Nsam,
			 char **var_names, char *file_name)
{
  char nameout[400]; int i, j;
  sprintf(nameout, "%s_matrix.dat", file_name);
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "# ");
  for(i=0; i<Nvar; i++)fprintf(file_out, " %s", var_names[i]);
  fprintf(file_out, "\n");
  for(j=0; j<Nsam; j++){
    for(i=0; i<Nvar; i++)fprintf(file_out, "%.3f ", data[i][j]);
    fprintf(file_out, "\n");
  }
  printf("Writing matrix with selected variables in %s\n", nameout);
  fclose(file_out);
}
 
int Select_points(float **data_select, float **data, int *index,
		  int *cluster, int clust, int Nsam, int Nvar)
{
  int i, j, n=0;
  for(i=0; i<Nsam; i++){
    if(cluster[i]==clust){
      index[n]=i;
      for(j=0; j<Nvar; j++)data_select[j][n]=data[j][i];
      n++;
    }
  }
  return(n);
}

int Correlation(double **C_matrix, float **data, int Nvar, int Nsam)
{

  double *average=malloc(Nsam*sizeof(double));
  double *deviation=malloc(Nsam*sizeof(double));
  int k1, k2, i;
  for(k1=0; k1<Nvar; k1++){
    double ave=0, dev=0;
    float *dat1=data[k1];
    for(i=0; i<Nsam; i++)ave+=dat1[i];
    ave/=Nsam; average[k1]=ave;

    for(k2=0; k2<= k1; k2++){
      double corr=0;
      float *dat2=data[k2];
      for(i=0; i<Nsam; i++)corr+=dat1[i]*dat2[i];
      corr=(corr/Nsam-(ave*average[k2]));
      C_matrix[k1][k2]=corr; C_matrix[k2][k1]=corr;
    }
  }
  free(average); free(deviation);
  return(0);
}

int Compute_PC(float **x_VarSam, float **PCA_load, float *eigen_values,
	       float **data, int num, int Nvar, int NPC)
{
  int i, j, k;
  // Compute weights
  //float *w=malloc(NPC*sizeof(float));
  //double wsum=0;
  //for(j=0; j<NPC; j++)wsum+=eigen_values[j];
  //for(j=0; j<NPC; j++)w[j]=eigen_values[j]/wsum;

  for(j=0; j<NPC; j++){
    float *x=x_VarSam[j];
    for(i=0; i<num; i++){
      double xx=0;
      for(k=0; k<Nvar; k++)xx+=PCA_load[k][j]*data[k][i];
      x[i]=xx; //*w
    }
  }
  return(0);
}

float AIC(float lik, int Npara, int Nsam){
  float AIC=-2*lik+2*Npara+(2*Npara*(Npara+1))/(float)(Nsam-Npara-1);
  return(AIC);
}

float BIC(float lik, int Npara, int Nsam){
  float BIC=-2*lik+Npara*log((float)(Nsam));
  return(BIC);
}

void Print_clusters(int **clust_store, float **P_store, int Nsam,
		    char *file_name, int ncl, int ncl_min, int ncl_max)
{
  if(ncl==1)return;
  char nameout[400]; int i;
  sprintf(nameout, "%s_clusters.txt", file_name);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing clusters in %s (1st col= index, 2nd col= P)\n", nameout);
  int n1=ncl_min, n; if(n1==1)n1=2;
  int n2=ncl_max; if(ncl==ncl_max)n2--;
  printf("Columns= %d %d-%d\n", ncl, n2, n1);
  for(i=0; i<Nsam; i++){
    fprintf(file_out, "%d\t%.3f",clust_store[ncl][i], P_store[ncl][i]);
    for(n=n2; n>=n1; n--){
      fprintf(file_out, "\t%d\t%.3f",clust_store[n][i], P_store[n][i]);
    }
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}

Remove_file_extension(char *nameout, char *ext, char *filename){
  char *s1=nameout, *s2=filename;
  while((*s2!='\0')&&(*s2!='.')){*s1=*s2; s1++; s2++;}
  *s1='\0';
  s1=ext;
  while(*s2!='\0'){if(*s1!='.')*s1=*s2; s1++; s2++;} *s1='\0';
}

Remove_dir(char *nameout, char *filename){
  char *s=filename, *start=s;
  while(*s!='\0'){if(*s=='/')start=s+1; s++;}
  strcpy(nameout, start);
}


void Cluster_statistics(int *numclus, float **ave, float **dev, int ncl,
			int *cluster, float **x_VarSam, int Nsam,
			int j0, int Nvar, float *lcorr)
{
  int i, j, k;
  for(k=0; k<ncl; k++){
    numclus[k]=0;
    for(j=j0; j<j0+Nvar; j++){
      ave[k][j]=0; dev[k][j]=0;
    }
  }
  for(i=0; i<Nsam; i++){
    k=cluster[i];
    if((k<0)||(k>=ncl)){
      printf("ERROR, sample %d has wrong cluster index %d\n", i, k);
      continue;
    }
    numclus[k]++;
    float *a=ave[k]+j0, *d=dev[k]+j0; 
    for(j=0; j<Nvar; j++){
      float x=x_VarSam[j][i];
      *a+=x; *d+=x*x; a++; d++;
    }
  }
  for(k=0; k<ncl; k++){
    int n=numclus[k]; if(n<2)continue;
    float *a=ave[k]+j0, *d=dev[k]+j0;
    for(j=0; j<Nvar; j++){
      (*a)/=n;
      (*d)=(*d)/n-(*a)*(*a);
      (*d)=SIGMA*sqrt((*d)*lcorr[j+j0]/(n-1));
      a++; d++;
    }
  }
}


int Print_profile(int *numclus, float **ave, float **dev, char **var_names,
		  int nclust, int Nsam, int Nvar, char *nameout, char *open)
{
  int i, j, k;
  char filename[400];
  sprintf(filename, "%s_cluster_profiles.dat", nameout);
  printf("Writing %s\n", filename);
  FILE *file_out=fopen(filename, open);
  fprintf(file_out, "# %d clusters\n", nclust);
  fprintf(file_out, "#");
  for(k=1; k<=nclust; k++)fprintf(file_out, "  Clust%d se ", k);
  fprintf(file_out, " Var\n");
  for(j=0; j<Nvar; j++){
    for(k=0; k<nclust; k++)
      fprintf(file_out, "  %6.3f %5.3f", ave[k][j], dev[k][j]);
    fprintf(file_out, "  %s\n", var_names[j]);
  }
  for(k=0; k<nclust; k++){
    float p=numclus[k]/(float)Nsam;
    fprintf(file_out, " %6d %4.0f", numclus[k], SIGMA*sqrt(numclus[k]*(1-p)));
  }
  fprintf(file_out, "  Elements\n");
  for(k=0; k<nclust; k++){
    float p=numclus[k]/(float)Nsam;
    fprintf(file_out, " %6.1f %4.1f", p*100, SIGMA*sqrt(p*(1-p)/Nsam)*100);
  }
  fprintf(file_out, "  Percentage\n");
  fclose(file_out);
}

Print_sim(float **Clust_sim, int nclust, int Nvar, char *nameout,
	  char *open, char *what)
{
  int i, j, k;
  char filename[400];
  sprintf(filename, "%s_cluster_similarity.dat", nameout);
  printf("Writing %s\n", filename);
  FILE *file_out=fopen(filename, open);
  if(strcmp(open,"w")==0)
    fprintf(file_out,"# %d variables\n", Nvar);
  fprintf(file_out, "###     %d clusters\n", nclust);
  fprintf(file_out, "###     %s\n", what);
  fprintf(file_out, "###     ");
  for(k=1; k<=nclust; k++)fprintf(file_out, " Clust%d", k);
  fprintf(file_out, "\n");
  for(k=0; k<nclust; k++){
    fprintf(file_out, "Clust%d ", k+1);
    if(k<9)fprintf(file_out, " ");
    for(i=0; i<7*k; i++)fprintf(file_out, " ");
    for(j=k; j<nclust; j++)
      fprintf(file_out, " %6.3f", Clust_sim[k][j]);
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}

Split_file(char *file_split, int nclust, int *cluster, int Nsam){
  int i, j, k, n=0;
  FILE *file_in=fopen(file_split, "r");
  if(file_in==NULL){
    printf("ERROR, file to be split %s does not exist\n", file_split);
    return(-1);
  }
  char string[1000];
  while(fgets(string, sizeof(string), file_in)!=0){
    if(string[0]=='#')continue; n++;
  }
  fclose(file_in);
  if(n!=Nsam){
    printf("ERROR, number of objects in file %s not equal to %d\n",
	   file_split, nclust);
    return(-2);
  }
  char filename[400], nameout[400], ext[100];
  Remove_file_extension(nameout, ext, file_split);
  FILE **file_out=malloc(nclust*sizeof(FILE *));

  for(k=0; k<nclust; k++){
    sprintf(filename, "%s%d%s", nameout, k+1, ext);
    printf("Writing %s\n", filename);
    file_out[k]=fopen(filename, "w");
  }
  file_in=fopen(file_split, "r"); n=0;
  while(fgets(string, sizeof(string), file_in)!=0){
    if(string[0]=='#')continue;
    k=cluster[n]; n++;
    if((k<0)||(k>=nclust)){
      printf("ERROR, element %d cluster %d out of bound\n", n-1, k+1);
      continue;
    }
    fprintf(file_out[k], "%s", string);
  }
  fclose(file_in);
  for(k=0; k<nclust; k++)fclose(file_out[k]);
}

int Print_cluster_PC(float **x_VarSam, int *cluster,
		     int nclust, int Nsam, int Npca, char *nameout)
{
  int i, j, k;
  char filename[400];
  sprintf(filename, "%s_cluster_PC.dat", nameout);
  printf("Writing %s\n", filename);
  FILE *file_out=fopen(filename, "w");
  fprintf(file_out, "# ");
  for(k=0; k<nclust; k++){
    fprintf(file_out, "# Clust%d\n", k);
    for(i=0; i<Nsam; i++){
      if(cluster[i]!=k)continue;
      for(j=0; j<Npca; j++)
	fprintf(file_out, "%.3f ",x_VarSam[j][i]);
      fprintf(file_out, "\n");
    }
    fprintf(file_out, "&\n");
  }
  fclose(file_out);
}

char **Set_names(char **var_names, int Nvar, int NPCA){
  char **names=malloc((Nvar+NPCA)*sizeof(char *)); int i;
  for(i=0; i<(Nvar+NPCA); i++){
    names[i]=malloc(80*sizeof(char));
    if(i<Nvar){strcpy(names[i], var_names[i]);}
    else{sprintf(names[i], "PC%d", i-Nvar+1);}
  }
  return(names);
}

static int Print_parameters(int *numcl,
			    float *tau, float **corr,
			    float **mu, float ***sig,
			    int ncl, int Nvar, FILE *file)
{
  int k, j, i;
  char txt[1000000];
  sprintf(txt, "");
  for(k=0; k<ncl; k++){
    sprintf(txt, "%sClus%d %4d samples, ", txt, k+1, numcl[k]);
    sprintf(txt, "%stau= %.2f mu= ", txt, tau[k]);
    for(j=0; j<Nvar; j++)sprintf(txt, "%s%6.3f ", txt, mu[k][j]);
    sprintf(txt, "%s sig= ", txt);
    for(j=0; j<Nvar; j++)sprintf(txt, "%s%5.3f ", txt, sig[k][j][j]);
    sprintf(txt, "%s corr= ", txt);
    for(j=0; j<ncl; j++)sprintf(txt, "%s%4.2f ", txt, corr[k][j]);
    sprintf(txt, "%s\n", txt);
  }
  if(file==NULL){
    printf("ERROR, can not print parameters in file:\n");
    printf("%s", txt);
  }else{
    fprintf(file, "%s", txt);
  }
}

static int Print_propensity(float *tau, float **corr, int ncl, FILE *file)
{
  char txt[1000000];
  sprintf(txt, "");
  int k, k1;
  sprintf(txt, "%s     ", txt);
  for(k=0; k<ncl; k++)sprintf(txt, "  %scl. %2d", txt, k);
  sprintf(txt, "%s\n", txt);

  for(k=0; k<ncl; k++){
    sprintf(txt, "%scl. %2d ", txt, k);
    for(k1=0; k1<ncl; k1++){
      sprintf(txt, " %.3f", txt, log(corr[k][k1]-log(tau[k1])));
    }
    sprintf(txt, "%s\n", txt);
  }
  
  if(file==NULL){
    printf("ERROR, can not print propensity in file\n");
    printf("%s", txt);
  }else{
    fprintf(file, "%s", txt);
  }

}

void Get_lcorr(float *lcorr, float **data, int N, int k0, int Nvar)
{
  int k, i, l, L=30; if(L>(N-1))L=N-1;
  float *x=malloc(L*sizeof(float));
  for(l=0; l<L; l++)x[l]=l;
  float *y=malloc(L*sizeof(float));
  for(k=0; k<Nvar; k++){
    double sum=0; int lmax=0;
    float *d=data[k], r, offset, slope;
    for(l=0; l<L; l++){
      double c=0; for(i=0; i<N-l; i++)c+=d[i]*d[i+l];
      c/=(N-l);
      if(c>0){y[l]=log(c);}else{break;}
    }
    lmax=l;
    r=Corr_coeff(&slope, &offset, x, y, lmax);
    slope = -slope;
    if((slope > 0)&&(slope < 1)){lcorr[k0+k]=1./slope;}
    else{lcorr[k0+k]=1;}
  }
}

float Corr_coeff(float *slope, float *offset, float *xx, float *yy, int n)
{
  float *x=xx, *y=yy; int i;
  double x1=0, x2=0, y1=0, y2=0, xy=0;
  for(i=0; i<n; i++){
    x1+=(*x); x2+=(*x)*(*x);
    y1+=(*y); y2+=(*y)*(*y);
    xy+=(*x)*(*y); x++; y++;
  }
  x1/=n; y1/=n;
  xy=(xy/n-x1*y1);
  x2=(x2/n-x1*x1);
  y2=(y2/n-y1*y1);
  *slope = xy/x2;
  *offset=(y1-(*slope)*x1);
  return(xy/sqrt(x2*y2));
}
