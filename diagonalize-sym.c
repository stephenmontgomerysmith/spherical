#include "spherical.h"

/* 
 * Find the eigenvalues of a symmetric matrix A using Jacobi Method.  The
 * eigenvalues are returned in eval[0],...,eval[n-1].  The eigenvectors
 * are returned in evec[0],...,evec[n-1], and are each norm one.  The
 * algorithm destroys the entries of A.
 *
 * The Jacobi Method works by repeatedly diagonalizing each 2 by 2 diagonal
 * sub-matrix.  Since orthogonal changes of basis preserve the sum of the
 * squares of the entries of the matrix, it is not hard to see that the sum of
 * the squares of the off diagonal entries converge to zero.  It is known to
 * converge extremely quickly.
 *
 * The main reasons I use a home made program instead of a ready made program
 * are these:
 * 1. lapack is in Fortran, and this is difficult (although not impossible)
 *    to link with C programs;
 * 2. gsl is GPL, and requires that this software be GPL if I use it;
 * 3. I don't see how to use either with mingw-32, and so I would be unable
 *    to cross compile for windows under FreeBSD;
 * 4. This algorithm doesn't need to be particularly fast, as computing
 *    eigenvalues and eigenvectors is not a big bottleneck for this suite of
 *    programs.
 */

void diagonalize_sym(int n, REAL *A, REAL *eval, REAL *evec) {
  int i, j, k, done;
  REAL discr, aii, aij, ajj, aik, ajk, v1, v2, norm;

  for (i=0;i<n;i++) for (j=0;j<n;j++)
    evec[i*n+j] = (i==j)?1:0;
  do {
    done = 1;
    for (i=0;i<n;i++) for (j=0;j<i;j++) {
/* Diagonalize the diagonal 2 by 2 sub-matrix along the ith and jth axes. */
      aii = A[i*n+i];
      aij = A[i*n+j];
      ajj = A[j*n+j];
      if (fabs(aij)>1e-100) {
        discr = hypot(aii-ajj,2*aij);
/* Compute eigenvalues of 2 by 2 sub-matrix. */
        if (aii+ajj>0) {
          A[i*n+i] = (aii+ajj+discr)/2;
          A[j*n+j] = (aii*ajj-aij*aij)/A[i*n+i];
        } else {
          A[j*n+j] = (aii+ajj-discr)/2;
          A[i*n+i] = (aii*ajj-aij*aij)/A[j*n+j];
        }
        A[i*n+j] = A[j*n+i] = 0;
/* Compute normalized eigenvector corresponding to first eigenvalue. */
        if (aii>ajj) {
          v1 = (aii-ajj+discr)/2;
          v2 = aij;
        } else {
          v1 = aij;
          v2 = (ajj-aii+discr)/2;
        }
        norm = hypot(v1,v2);
        v1 /= norm;
        v2 /= norm;
/* Apply change of basis to the rest of the matrix. */
        for (k=0;k<n;k++) if (k!=i && k!=j) {
          aik = A[i*n+k];
          ajk = A[j*n+k];
          A[i*n+k] = A[k*n+i] = v1*aik+v2*ajk;
          A[j*n+k] = A[k*n+j] = -v2*aik+v1*ajk;
        }
/* Apply change of basis to evec. */
        for (k=0;k<n;k++) {
          aik = evec[i*n+k];
          ajk = evec[j*n+k];
          evec[i*n+k] = v1*aik+v2*ajk;
          evec[j*n+k] = -v2*aik+v1*ajk;
        }
        done = 0;
      }
    }
  } while (!done);
  for (i=0;i<n;i++)
    eval[i] = A[i*n+i];
}

/* If you want to test it:

#include <stdio.h>

main () {
  REAL A[9] = {2,3,4,3,4,5,4,5,6}};
  REAL evec[9];
  int i,j,k;
  REAL s;

  diagonalize_sym(3,A,eval,evec);

  for (j=0;j<3;j++) {
    for (k=0;k<3;k++) {
      s = 0;
      for (i=0;i<3;i++)
        s += eval[i]*evec[i*3+j]*evec[i*3+k];
      printf("%g ",s);
    }
    printf("\n");
  }
}

*/
