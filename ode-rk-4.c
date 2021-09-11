#include "spherical.h"

static REAL *k1, *k2, *k3, *k4;
static REAL *tempx;
static int first = 1;

void ode_rk_4_solve(REAL *t, REAL *x, REAL h,
                    void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                    param_list_t *param) {
  int i;

  if (first) {
    first = 0;
    k1 = (REAL*)malloc(sizeof(REAL)*param->length);
    k2 = (REAL*)malloc(sizeof(REAL)*param->length);
    k3 = (REAL*)malloc(sizeof(REAL)*param->length);
    k4 = (REAL*)malloc(sizeof(REAL)*param->length);
    tempx = (REAL*)malloc(sizeof(REAL)*param->length);
  }
  derivs(*t,x,k1,param);
  for (i=0;i<param->length;i++) tempx[i] = x[i] + h*k1[i]/2;
  derivs(*t+h/2,tempx,k2,param);
  for (i=0;i<param->length;i++) tempx[i] = x[i] + h*k2[i]/2;
  derivs(*t+h/2,tempx,k3,param);
  for (i=0;i<param->length;i++) tempx[i] = x[i] + h*k3[i];
  derivs(*t+h,tempx,k4,param);
  for (i=0;i<param->length;i++) x[i] += h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
  *t += h;
}
