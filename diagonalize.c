#include <malloc.h>

/* nrutil.h/nrutil.c */
extern double *dvector( long, long);
extern double **dmatrix( long, long, long, long);
extern void free_dmatrix( double **, long, long, long, long);
extern void free_dvector( double *, long, long);
/* tred2.c (float -> double) */
extern void tred2( double **, int, double *, double *);
/* tqli.c (float -> double, bolean return value)  */
extern int tqli( double *, double *, int, double **);
/* jacobi.c (float -> double)  */
extern void jacobi( double **, int, double *, double **, int *);
/* eigsrt.c (float -> double)  */
extern void eigsrt( double *, double **, int);

void
Diagonalize(int N, double **MATRIX, float *eigen_values, float **eigen_vector)
{
  double **a_matrix, **v_matrix, **eigen_vt, *eigen_va, *e_vector;
  int jacobi_rotations=0, i, j;

  e_vector= dvector( (long)1, (long)N);
  eigen_va= dvector( (long)1, (long)N);
  a_matrix= dmatrix( (long)1, (long)N, (long)1, (long)N);
  v_matrix= dmatrix( (long)1, (long)N, (long)1, (long)N);

  /* fill matrix a_matrix */
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      a_matrix[i+1][j+1]=MATRIX[i][j];
    }
  }

  /* obtain eigensystem */
  tred2( a_matrix, (int)N, eigen_va, e_vector);
  if( tqli(eigen_va, e_vector, (int)N, a_matrix))
    {
      eigen_vt=a_matrix;
    }
  else
    {
      /* tqli failed */
      /* reconstruct matrix a_matrix (destroyed by tred2/tqli!) */
      for(i=0; i<N; i++){
	for(j=0; j<N; j++){
	  a_matrix[i+1][j+1]=MATRIX[i][j];
	}
      }

      /* obtain eigensystem */
      jacobi( a_matrix, (int)N, eigen_va, v_matrix, &jacobi_rotations);
      eigen_vt=v_matrix;
    }

  eigsrt(eigen_va, eigen_vt, (int)N);
  for(i=0; i<N; i++){
    eigen_values[i]=eigen_va[i+1];
    for(j=0; j<N; j++){
      eigen_vector[i][j]=eigen_vt[i+1][j+1];
    }
  }

  free_dvector(e_vector, 1, N);
  free_dvector(eigen_va, 1, N);
  free_dmatrix(a_matrix, 1, N, 1, N);
  free_dmatrix(v_matrix, 1, N, 1, N);

  return;
}


void
f_Diagonalize(int N, float **MATRIX, float *eigen_values, float **eigen_vector)
{
  double **a_matrix, **v_matrix, **eigen_vt, *eigen_va, *e_vector;
  int jacobi_rotations=0, i, j;

  e_vector= dvector( (long)1, (long)N);
  eigen_va= dvector( (long)1, (long)N);
  a_matrix= dmatrix( (long)1, (long)N, (long)1, (long)N);
  v_matrix= dmatrix( (long)1, (long)N, (long)1, (long)N);

  /* fill matrix a_matrix */
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      a_matrix[i+1][j+1]=MATRIX[i][j];
    }
  }

  /* obtain eigensystem */
  tred2( a_matrix, (int)N, eigen_va, e_vector);
  if( tqli(eigen_va, e_vector, (int)N, a_matrix))
    {
      eigen_vt=a_matrix;
    }
  else
    {
      /* tqli failed */
      /* reconstruct matrix a_matrix (destroyed by tred2/tqli!) */
      for(i=0; i<N; i++){
	for(j=0; j<N; j++){
	  a_matrix[i+1][j+1]=MATRIX[i][j];
	}
      }

      /* obtain eigensystem */
      jacobi( a_matrix, (int)N, eigen_va, v_matrix, &jacobi_rotations);
      eigen_vt=v_matrix;
    }

  eigsrt(eigen_va, eigen_vt, (int)N);
  for(i=0; i<N; i++){
    eigen_values[i]=eigen_va[i+1];
    for(j=0; j<N; j++){
      eigen_vector[i][j]=eigen_vt[i+1][j+1];
    }
  }

  free_dvector(e_vector, 1, N);
  free_dvector(eigen_va, 1, N);
  free_dmatrix(a_matrix, 1, N, 1, N);
  free_dmatrix(v_matrix, 1, N, 1, N);

  return;
}


