#include "spherical.h"

static REAL *k1, *k2, *k3, *k4;
static REAL *diffx, *olddiffx1, *olddiffx2, *olddiffx3;
static REAL *tempx;
static int iter = 0;

void ode_adams_bash_4_solve(REAL *t, REAL *x, REAL h,
                            void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                            param_list_t *param) {
  int i;

  if (iter==0) {
    k1 = (REAL*)malloc(sizeof(REAL)*param->length);
    k2 = (REAL*)malloc(sizeof(REAL)*param->length);
    k3 = (REAL*)malloc(sizeof(REAL)*param->length);
    k4 = (REAL*)malloc(sizeof(REAL)*param->length);
    tempx = (REAL*)malloc(sizeof(REAL)*param->length);
    diffx = (REAL*)malloc(sizeof(REAL)*param->length);
    olddiffx1 = (REAL*)malloc(sizeof(REAL)*param->length);
    olddiffx2 = (REAL*)malloc(sizeof(REAL)*param->length);
    olddiffx3 = (REAL*)malloc(sizeof(REAL)*param->length);
  }
  if (iter<=2) {
/* Runge-Kutta method of order 4. */
    memcpy(olddiffx3,olddiffx2,sizeof(REAL)*param->length);
    memcpy(olddiffx2,olddiffx1,sizeof(REAL)*param->length);
    derivs(*t,x,k1,param);
    memcpy(olddiffx1,k1,sizeof(REAL)*param->length);
    for (i=0;i<param->length;i++) tempx[i] = x[i] + h*k1[i]/2;
    derivs(*t+h/2,tempx,k2,param);
    for (i=0;i<param->length;i++) tempx[i] = x[i] + h*k2[i]/2;
    derivs(*t+h/2,tempx,k3,param);
    for (i=0;i<param->length;i++) tempx[i] = x[i] + h*k3[i];
    derivs(*t+h,tempx,k4,param);
    for (i=0;i<param->length;i++)
      x[i] += h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
  } else {
/* Adams-Bashforth method of order 4. */
    derivs(*t,x,diffx,param);
    for (i=0;i<param->length;i++)
      x[i] += h*(55.0/24*diffx[i]-59.0/24*olddiffx1[i]+37.0/24*olddiffx2[i]-3.0/8*olddiffx3[i]);
    memcpy(olddiffx3,olddiffx2,sizeof(REAL)*param->length);
    memcpy(olddiffx2,olddiffx1,sizeof(REAL)*param->length);
    memcpy(olddiffx1,diffx,sizeof(REAL)*param->length);
/* Use this instead of the three lines above if you want to exchange pointers
   to data instead of data itself.  (But it didn't seem to make a big speed
   difference.)
    REAL *temp;
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
