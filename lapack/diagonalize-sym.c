#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Interface to lapack Fortran routine. */

static double *work;
static int ldwork;
static int last_n = 0;

void dsyev_(char *jobz, char *uplo, int *n, double *A, int *lda,
            double *w, double *work, int *ldwork, int *info);

void diagonalize_sym(int n, double A[n][n], double eval[n], double evec[n][n]) {
  char jobz='V', uplo='U';
  int info;
  double ret_ldwork;

  if (n!=last_n) {
    ldwork = -1;
    dsyev_(&jobz, &uplo, &n, (double*)A, &n, eval, &ret_ldwork, &ldwork, &info);
    if (info != 0) {
      fprintf(stderr,"Problem computing eigenvalues.\n");
      exit(1);
    }
    ldwork = (int)ret_ldwork;
    if (last_n==0)
      work = malloc(sizeof(double)*ldwork);
    else
      work = realloc(work,sizeof(double)*ldwork);
    last_n = n;
  }

  dsyev_(&jobz, &uplo, &n, (double*)A, &n, eval, work, &ldwork, &info);
  if (info != 0) {
    fprintf(stderr,"Problem computing eigenvalues.\n");
    exit(1);
  }
  memcpy(evec,A,sizeof(double[n][n]));
}

/* If you want to test it:

#include <stdio.h>

int main () {
  double A[3][3] = {{2,13,4},
                    {13,4,55},
                    {4,55,86}};
  double evec[3][3], eval[3];
  int i,j,k;
  double s;

  diagonalize_sym(3,A,eval,evec);

  for (j=0;j<3;j++) {
    for (k=0;k<3;k++) {
      s = 0;
      for (i=0;i<3;i++)
        s += eval[i]*evec[i][j]*evec[i][k];
      printf("%g ",s);
    }
    printf("\n");
  }
  return 0;
}

*/
