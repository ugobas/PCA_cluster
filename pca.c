#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#define NS_DEF 0

/* Needs name of input file
   Format of the input file: ### with the names of the n variables
   L lines with n variables per line. */

void help(char **argv);
void Get_args(char *file_name, char *name_out, int *N_SHOW, int *Print,
	      int *MATRIX, int argc, char **argv);

extern void
Diagonalize(int N, double **MATRIX,
	    float *eigen_values, float **eigen_vector);

float **Eigensystem(int, float **, double **);
float **Read_data(char *file_name, int *N, int *N_data, char ***, char ***,
		  int **select_var);
double **Compute_corr(int N, int N_data, float **data, float **,
		      double **, double **);
int Print_results(float *, float **, float **, double **,
		  char **, char **, int, int, char *, double *, double *,
		  int *N_SHOW, int Print);
int Read_var(char *string, int N, int *variables, char **var_names);


main(int argc, char **argv){
  int N_SHOW=NS_DEF, Print=0;
  float **data, **data_norm, **eigen_vector=NULL, *eigen_values=NULL;
  char **var_names, **ind_names, file_name[200], name_out[200], name[200];
  double *average, *deviation;
  double **correlation;
  int N=0, N_data, i, j;
  int *select_var;
  int MATRIX;

  Get_args(file_name, name_out, &N_SHOW, &Print, &MATRIX, argc, argv);
  data=Read_data(file_name, &N, &N_data, &var_names, &ind_names, &select_var);
  N=Select_data(data, N, N_data, var_names, select_var);

  data_norm=malloc(N*sizeof(float *));
  for(i=0; i<N; i++)data_norm[i]=malloc(N_data*sizeof(float));
  correlation=Compute_corr(N, N_data, data, data_norm, &average, &deviation);
  printf("%d elements, %d selected variables:\n", N_data, N);
  float norm=(3*2.0)/sqrt(N_data);
  for(i=0; i<N; i++){
    printf("%9.3g\t%.1g\t%s\n", average[i], norm*deviation[i], var_names[i]);
  }

  eigen_vector=Eigensystem(N, &eigen_values, correlation);
  Print_results(eigen_values, eigen_vector, data_norm, correlation,
		var_names, ind_names, N, N_data, name_out,
		average, deviation, &N_SHOW, Print);

  if(MATRIX){
    char nameout[400];
    sprintf(nameout, "%s_matrix.dat", file_name);
    FILE *file_out=fopen(nameout, "w");
    fprintf(file_out, "# ");
    for(i=0; i<N; i++)fprintf(file_out, " %s", var_names[i]);
    fprintf(file_out, "\n");
    for(j=0; j<N_data; j++){
      for(i=0; i<N; i++)fprintf(file_out, "%.3f ", data[i][j]);
      fprintf(file_out, "\n");
    }
    printf("Writing matrix with selected variables in %s\n", nameout);
    fclose(file_out);
  }

  return(0);
}

void Get_args(char *file_name, char *name_out, int *N_SHOW, int *Print,
	      int *MATRIX, int argc, char **argv)
{
  int i;
  if(argc<2)help(argv);
  for(i=1; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0){
      help(argv);
    }else if(i==1){
      strcpy(file_name, argv[1]);
      sprintf(name_out, "%s.pca\0", file_name);
    }else if(strncmp(argv[i], "-n", 2)==0){
      i++; sscanf(argv[i], "%d", N_SHOW);
      if(*N_SHOW >0){
	printf("Number of modes to print: %d\n", *N_SHOW);
      }else{
	printf("WARNING, %d not allowed as number of modes to print.",*N_SHOW);
	printf("set to default %d\n", NS_DEF);
	*N_SHOW=NS_DEF;
      }
    }else if(strncmp(argv[i], "-print", 6)==0){
      *Print=1;
    }else if(strncmp(argv[i], "-matrix", 6)==0){
      *MATRIX=1;
    }else{
      printf("WARNING, unknown option %s\n", argv[i]);
    }
  }
  if(file_name[0]==0)help(argv);
}

void help(char **argv){
  printf("Program %s, author Ugo Bastolla.\n", argv[0]);
  printf("Usage: %s <input_file> (mandatory)\n", argv[0]);
  printf("Options:\n");
  printf("   -h help\n");
  printf("   -n <PCA loads to print>\n");
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
  *N=0; *N_data=0;

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
    sprintf((*var_names)[i], "X%d\0", i); 
  }
  *ind_names=malloc((*N_data)*sizeof(char *));
  for(i=0; i< (*N_data); i++){
    (*ind_names)[i]=malloc(100*sizeof(char));
    sprintf((*ind_names)[i], "I%d\0", i);
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
	//string[0]=' '; sscanf(string, "%s", (*ind_names)[m_data-1]);
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

int
Print_results(float *eigen_values, float **eigen_vector, float **data_norm,
	      double **correlation, char **var_names, char **ind_names,
	      int N, int N_data, char *name_out, double *average,
	      double *deviation, int *N_SHOW, int Print)
{
  if(*N_SHOW>=N)*N_SHOW=N-1;
  FILE *file_out=fopen(name_out, "w");
  int v1, v2, e1, i;  double tot=0, e;
  float r1, r2, *mean=malloc(N*sizeof(float));

  printf("Writing %s\n", name_out);
  fprintf(file_out, "# %d variables, %d elements %d PC\n",
	  N, N_data, *N_SHOW);

  fprintf(file_out, "# Correlation matrix\n");
  for(v1=0; v1<N; v1++){
    /*for(v2=v1+1; v2<N; v2++){
      fprintf(file_out, "# %6.3f ", correlation[v1][v2]);
      fprintf(file_out, "%12s %12s\n", var_names[v1], var_names[v2]);
      }*/
    for(v2=v1+1; v2<N; v2++)
      fprintf(file_out, "%s\t%s\t%.3f\n",
	      var_names[v1], var_names[v2], correlation[v1][v2]);
  }
  if(eigen_values==NULL){fclose(file_out); return(0);}

  int N_sig=0;
  fprintf(file_out, "# Spectrum: n lambda_n cumulative\n");
  for(e1=0; e1<N; e1++){
    tot+=eigen_values[e1];
    e=0; for(i=0; i<N; i++)e+=eigen_vector[i][e1]; e/=N;
    if(e<0){
      for(i=0; i<N; i++)eigen_vector[i][e1]=-eigen_vector[i][e1]; e=-e;
    }
    mean[e1]=e;
    fprintf(file_out, "# %2d %6.3f %.3f\n", e1, eigen_values[e1], tot/N);
    if(eigen_values[e1]>1.0)N_sig++;
  }
  if(*N_SHOW==0)*N_SHOW=N_sig;

  fprintf(file_out, "# %d variables, %d elements\n", N, N_data);
  fprintf(file_out, "#      Ave       s.e.       sigma  ");
  for(e1=0; e1<*N_SHOW; e1++)fprintf(file_out, " PC%d   ", e1);
  fprintf(file_out, "sum(PC)^2  var.\n");
  for(v1=0; v1<N; v1++){
    fprintf(file_out, "# %9.3g %9.3g %9.3g   ",
	    average[v1], deviation[v1]/sqrt(N), deviation[v1]);
    tot=0;
    for(e1=0; e1<*N_SHOW; e1++){
      if(e1>=N)break;
      r1= eigen_vector[v1][e1]*sqrt(eigen_values[e1]); tot+=r1*r1;
      fprintf(file_out, "%6.3f ", r1);
    }
    fprintf(file_out, " %.3f %s\n", tot, var_names[v1]);
  }

  if(Print==0)return(0);
  fprintf(file_out, "# Components 1...%d  Individual\n", *N_SHOW);
  for(i=0; i<N_data; i++){
    for(e1=0; e1<*N_SHOW; e1++){     
      if(e1>=N)break; r1=0;
      for(v1=0; v1<N; v1++)
	r1+=eigen_vector[v1][e1]*(data_norm[v1][i]);
      fprintf(file_out, "%6.3f ", r1);
    }
    fprintf(file_out, "  %s\n", ind_names[i]);
  }

  return(0);
}

double **
Compute_corr(int N, int N_data, float **data, float **data_norm,
	     double **average, double **deviation)
{
  double **C_matrix=malloc(N*sizeof(double *));
  double corr=0, ave, dev;
  float *dat1, *dat2;
  int k1, k2, i;

  for(k1=0; k1<N; k1++){
    C_matrix[k1]=malloc(N*sizeof(double));
  }

  *average=malloc(N*sizeof(double));
  *deviation=malloc(N*sizeof(double));
  for(k1=0; k1<N; k1++){
    ave=0; dev=0;
    dat1=data[k1]; C_matrix[k1][k1]=1.;
    for(i=0; i<N_data; i++){
      ave+=dat1[i]; dev+=dat1[i]*dat1[i];
    }   
    ave/=N_data; dev=(dev/N_data)-(ave*ave); dev=sqrt(dev);
    for(i=0; i<N_data; i++){
      data_norm[k1][i]=(dat1[i]-ave)/dev;
    }
    (*average)[k1]=ave; (*deviation)[k1]=dev;

    for(k2=0; k2< k1; k2++){
      corr=0; dat2=data[k2];
      for(i=0; i<N_data; i++){
	corr+=dat1[i]*dat2[i];
      }
      corr=(corr/N_data-(ave*(*average)[k2]))/
	(dev*(*deviation)[k2]);

      C_matrix[k1][k2]=corr; C_matrix[k2][k1]=corr;
    }
  }

  return(C_matrix);
}

float **
Eigensystem(int N, float **eigen_values, double **matrix)
{
  int i;
  float **eigen_vector=malloc(N*sizeof(float *));

  (*eigen_values)=malloc(N*sizeof(float));
  for(i=0; i<N; i++)eigen_vector[i]=malloc(N*sizeof(float));
  Diagonalize(N, matrix, *eigen_values, eigen_vector);
  return(eigen_vector);
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
