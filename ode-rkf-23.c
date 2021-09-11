#include "spherical.h"

/* Algorithm from "Numerical Recipes in C" by Press, Teukolsky, Vetterling, and Flannery;
   numbers taken from ode23.m in odepkg package from octave forge. */

static REAL **k;
static REAL *tempx;
static int first = 1;

static REAL
a[3] = {0, 1./2, 1},
b[3][3] = {{0, 0, 0},
{1./2, 0, 0},
{-1, 2, 0}},
c3[3] = {1./6, 2./3, 1./6},
/* c2[3] = {0, 1, 0}, */
c23[3] = {1./6, -1./3, 1./6}; /* c3-c2 */

void ode_rkf_23_solve(REAL *t, REAL *x, REAL *h_use,
                      void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                      param_list_t *param) {
  int i,j,m;
  REAL error, h, delta, temp1;

  if (first) {
    first = 0;
    k = (REAL**)malloc(sizeof(REAL*)*3);
    for (i=0;i<3;i++)
      k[i] = (REAL*)malloc(sizeof(REAL)*param->length);
    tempx = (REAL*)malloc(sizeof(REAL)*param->length);
  }

  h = *h_use;
  while (1) {
    for (i=0;i<3;i++) {
      memcpy(tempx,x,sizeof(REAL)*param->length);
      for (j=0;j<i;j++)
        for (m=0;m<param->length;m++)
          tempx[m] += h*b[i][j]*k[j][m];
      derivs(*t+a[i]*h,tempx,k[i],param);
    }
/* l_infinity error between RK2 and RK3 */
    error = 0;
    for (m=0;m<param->length;m++) {
      temp1 = 0;
      for (i=0;i<3;i++)
        temp1 += c23[i]*k[i][m];
      if (fabs(temp1)>error) error = fabs(temp1);
    }
    error = h*error/param->tol;
    if (error<=1) break;
    delta = 0.9*pow(error,-0.5);
    if (delta>0.1)
      h = delta*h;
    else
      h = 0.1*h;
  }

/* RK order 3 */
  for (i=0;i<3;i++)
    for (m=0;m<param->length;m++)
      x[m] += h*c3[i]*k[i][m];
  *t += h;

  if (error>5.832e-3) /* (5/.9)^-3 */
    *h_use = 0.9*h*pow(error,-1./3);
  else
    *h_use = 5*h;
}
