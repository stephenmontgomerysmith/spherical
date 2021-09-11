#include "spherical.h"

static COMPLEX *diffx;
static COMPLEX *olddiffx;
static COMPLEX *tempx;
static int first = 1;

extern int max_order;

void ode_adams_bash_2_solve(double *t, COMPLEX *x, double h,
                            void derivs(double t, COMPLEX *x, COMPLEX *diffx)) {
  int i;

  if (first) {
    first = 0;
    diffx = malloc(sizeof(COMPLEX)*length);
    tempx = malloc(sizeof(COMPLEX)*length);
    olddiffx = malloc(sizeof(COMPLEX)*length);
/* Midpoint method */
    derivs(*t,x,diffx);
    memcpy(olddiffx,diffx,length*sizeof(COMPLEX));
    for (i=0;i<length;i++) tempx[i] = x[i] + h*diffx[i]/2;
    derivs(*t+h/2,tempx,diffx);
    for (i=0;i<length;i++) x[i] += h*diffx[i];
  } else {
/* Adams-Bashforth method of order 2 */
    derivs(*t,x,diffx);
    for (i=0;i<length;i++)
      x[i] = x[i] + h/2 * (3*diffx[i] - olddiffx[i]);
    memcpy(olddiffx,diffx,length*sizeof(COMPLEX));
  }
  *t += h;
}
