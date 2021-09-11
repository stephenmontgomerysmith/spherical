#include "spherical.h"

static COMPLEX *k1, *k2, *k3, *k4;
static COMPLEX *tempx;
static int first = 1;

extern int max_order;

void ode_rk_4_solve(double *t, COMPLEX *x, double h,
                    void derivs(double t, COMPLEX *x, COMPLEX *diffx)) {
  int i;

  if (first) {
    first = 0;
    k1 = malloc(sizeof(COMPLEX)*length);
    k2 = malloc(sizeof(COMPLEX)*length);
    k3 = malloc(sizeof(COMPLEX)*length);
    k4 = malloc(sizeof(COMPLEX)*length);
    tempx = malloc(sizeof(COMPLEX)*length);
  }
  derivs(*t,x,k1);
  for (i=0;i<length;i++) tempx[i] = x[i] + h*k1[i]/2;
  derivs(*t+h/2,tempx,k2);
  for (i=0;i<length;i++) tempx[i] = x[i] + h*k2[i]/2;
  derivs(*t+h/2,tempx,k3);
  for (i=0;i<length;i++) tempx[i] = x[i] + h*k3[i];
  derivs(*t+h,tempx,k4);
  for (i=0;i<length;i++) x[i] += h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
  *t += h;
}
