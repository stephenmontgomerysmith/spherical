#include "spherical.h"

static COMPLEX *k1, *k2, *k3, *k4;
static COMPLEX *diffx, *olddiffx1, *olddiffx2, *olddiffx3;
static COMPLEX *tempx;
static int iter = 0;

extern int max_order;

void ode_adams_bash_4_solve(double *t, COMPLEX *x, double h,
                            void derivs(double t, COMPLEX *x, COMPLEX *diffx)) {
  int i;

  if (iter==0) {
    k1 = malloc(sizeof(COMPLEX)*length);
    k2 = malloc(sizeof(COMPLEX)*length);
    k3 = malloc(sizeof(COMPLEX)*length);
    k4 = malloc(sizeof(COMPLEX)*length);
    tempx = malloc(sizeof(COMPLEX)*length);
    diffx = malloc(sizeof(COMPLEX)*length);
    olddiffx1 = malloc(sizeof(COMPLEX)*length);
    olddiffx2 = malloc(sizeof(COMPLEX)*length);
    olddiffx3 = malloc(sizeof(COMPLEX)*length);
  }
  if (iter<=2) {
/* Runge-Kutta method of order 4. */
    memcpy(olddiffx3,olddiffx2,sizeof(COMPLEX)*length);
    memcpy(olddiffx2,olddiffx1,sizeof(COMPLEX)*length);
    derivs(*t,x,k1);
    memcpy(olddiffx1,k1,sizeof(COMPLEX)*length);
    for (i=0;i<length;i++) tempx[i] = x[i] + h*k1[i]/2;
    derivs(*t+h/2,tempx,k2);
    for (i=0;i<length;i++) tempx[i] = x[i] + h*k2[i]/2;
    derivs(*t+h/2,tempx,k3);
    for (i=0;i<length;i++) tempx[i] = x[i] + h*k3[i];
    derivs(*t+h,tempx,k4);
    for (i=0;i<length;i++)
      x[i] += h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
  } else {
/* Adams-Bashforth method of order 4. */
    derivs(*t,x,diffx);
    for (i=0;i<length;i++)
      x[i] += h*(55.0/24*diffx[i]-59.0/24*olddiffx1[i]+37.0/24*olddiffx2[i]-3.0/8*olddiffx3[i]);
    memcpy(olddiffx3,olddiffx2,sizeof(COMPLEX)*length);
    memcpy(olddiffx2,olddiffx1,sizeof(COMPLEX)*length);
    memcpy(olddiffx1,diffx,sizeof(COMPLEX)*length);
/* Use this instead of the three lines above if you want to exchange pointers
   to data instead of data itself.  (But it didn't seem to make a big speed
   difference.)
    COMPLEX *temp;
    temp = olddiffx3;
    olddiffx3 = olddiffx2;
    olddiffx2 = olddiffx1;
    olddiffx1 = diffx;
    diffx = temp;
*/
  }
  *t += h;
  iter++;
}
