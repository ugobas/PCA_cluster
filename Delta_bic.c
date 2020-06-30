#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "allocate.h"
#include "read_data.h"

int Get_para(char *file_name, int *NPC, int *Nsam);

int main(int argc, char **argv)
{
  if(argc<2){
    printf("ERROR, input file name must be specified\n");
  }
  char file_name[200];
  strcpy(file_name, argv[1]);

  int Nvar=0, Ndata=0, *select_var;
  char **var_names, **ind_names;
  float **data=Read_data(file_name, &Nvar, &Ndata,
			 &var_names, &ind_names, &select_var);

  char nameout[200], ext[80];
  Remove_file_extension(nameout, ext, file_name);
  sprintf(nameout, "%s_Delta_bic.dat", nameout);
  FILE *file_out=fopen(nameout, "w");
  printf("Writing %s\n", nameout);

  int NPC=0, Nsam;
  Get_para(file_name, &NPC, &Nsam);

  FILE *file_in=fopen(file_name, "r");
  char string[200];
  while(fgets(string, sizeof(string), file_in)!=0){
    if(string[0]=='#')fprintf(file_out, "%s", string);
  }
  fclose(file_in);

  int bic=3, imin=0, i, j;
  float min=data[bic][imin];
  float **y=Allocate_mat2_f(Ndata, Nvar);
  for(i=0; i<Ndata; i++){
    fprintf(file_out, "%2.0f", data[0][i]);
    for(j=1; j<Nvar; j++){
      y[i][j]=data[j][i]-data[j][0];
      //if(data[j][0])y[i][j]/=fabs(data[j][0]);
      if(j<4)y[i][j]/=NPC;
      fprintf(file_out, "\t%.4f", y[i][j]);
    }
    fprintf(file_out, "\n"); fflush(file_out);
    if(data[bic][i] < min){min=data[bic][i]; imin=i;}
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

int Get_para(char *file_name, int *NPC, int *Nsam){
  FILE *file_in=fopen(file_name, "r");
  char string[200]; int k=0;
  // # PCA_cluster Properties_genome_w360_s050_n4 4 PC 165478 samples
  char dum0[80], dum1[80], dum2[80], dum3[80]; 
  while(fgets(string, sizeof(string), file_in)!=0){
    if(string[0]=='#')k++;
    if(k==2){
      sscanf(string, "%s %s %s %d %s %d",
	     dum0, dum1, dum2, NPC, dum3, Nsam);
      break;
    }
  }
  fclose(file_in);
  printf("NPC= %d Nsam= %d\n", *NPC, *Nsam);
  return(0);
}
