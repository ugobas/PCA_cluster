#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "allocate.h"
#include "read_data.h"

int Nvar=16;

int Get_para(char *file_name, int *NPCA, int *Nsam);
void help(char **argv);
float AIC(float lik, int Npara, int Nsam);
float BIC(float lik, int Npara, int Nsam);

int main(int argc, char **argv)
{
  char file_name[200];
  if(argc<2){
    printf("ERROR, input file name must be specified\n");
    help(argv);
  }
  strcpy(file_name, argv[1]);

  
  int Nvar=0, N=0, *select_var;
  char **var_names, **ind_names;
  float **data=Read_data(file_name, &Nvar, &N,
			 &var_names, &ind_names, &select_var);

  // Corrections
  int NPCA=0, Nsam;
  Get_para(file_name, &NPCA, &Nsam);
  int i;
  for(i=0; i<N; i++){
    float lik =data[1][i]*NPCA*NPCA*Nsam;
    int N_para=data[0][i]*(NPCA + NPCA*(NPCA+1)/2-1)-1;
    data[1][i]=lik/Nsam; // likelihood
    data[2][i]=AIC(lik, N_para, Nsam)/Nsam;
    data[3][i]=BIC(lik, N_para, Nsam)/Nsam;
  }

  char nameout[200], ext[80];
  Remove_file_extension(nameout, ext, file_name);
  sprintf(nameout, "%s_Delta_bic.dat", nameout);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);



  FILE *file_in=fopen(file_name, "r");
  char string[200];
  while(fgets(string, sizeof(string), file_in)!=0){
    if(string[0]=='#')fprintf(file_out, "%s", string);
  }
  fclose(file_in);

  int bic=3, imin=0, j;
  float min=data[bic][imin];
  float **y=Allocate_mat2_f(N, Nvar);
  for(i=0; i<N; i++){
    fprintf(file_out, "%2.0f", data[0][i]);
    for(j=1; j<Nvar; j++){
      y[i][j]=data[j][i];
      y[i][j]-=data[j][0];
      //if(data[j][0])y[i][j]/=fabs(data[j][0]);
      fprintf(file_out, "\t%.4f", y[i][j]);
    }
    fprintf(file_out, "\n");
    if(y[i][bic] < min){min=y[i][bic]; imin=i;}
  }
  fprintf(file_out, "# Optimal BIC: %.0f %.4f\n",
	  data[0][imin], data[bic][imin]-data[bic][0]);
  fclose(file_out);

  printf("# Optimal BIC:\n%.0f", data[0][imin]);
  for(j=1; j<Nvar; j++)printf("\t%.4f", y[imin][j]);
  printf("\n");

  Empty_matrix_f(y, Nvar);
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

int Get_para(char *file_name, int *NPCA, int *Nsam){
  FILE *file_in=fopen(file_name, "r");
  char string[200]; int k=0;
  // # PCA_cluster Properties_genome_w360_s050_n4 4 PC 165478 samples
  char dum0[80], dum1[80], dum2[80], dum3[80]; 
  while(fgets(string, sizeof(string), file_in)!=0){
    if(string[0]=='#')k++;
    if(k==2){
      sscanf(string, "%s %s %s %d %s %d",
	     dum0, dum1, dum2, NPCA, dum3, Nsam);
      break;
    }
  }
  fclose(file_in);
  printf("NPCA= %d Nsam= %d\n", *NPCA, *Nsam);
  return(0);
}

void help(char **argv){
  printf("USAGE: %s <input_file>_bic.dat -n NPCA\n", argv[0]);
  exit(8);
}
