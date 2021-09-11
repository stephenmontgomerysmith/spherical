#include "spherical.h"

static REAL *diffx;
static REAL *olddiffx;
static REAL *tempx;
static int first = 1;

void ode_adams_bash_2_solve(REAL *t, REAL *x, REAL h,
                            void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                            param_list_t *param) {
  int i;

  if (first) {
    first = 0;
    diffx = (REAL*)malloc(sizeof(REAL)*param->length);
    tempx = (REAL*)malloc(sizeof(REAL)*param->length);
    olddiffx = (REAL*)malloc(sizeof(REAL)*param->length);
/* Midpoint method */
    derivs(*t,x,diffx,param);
    memcpy(olddiffx,diffx,param->length*sizeof(REAL));
    for (i=0;i<param->length;i++) tempx[i] = x[i] + h*diffx[i]/2;
    derivs(*t+h/2,tempx,diffx,param);
    for (i=0;i<param->length;i++) x[i] += h*diffx[i];
  } else {
/* Adams-Bashforth method of order 2 */
    derivs(*t,x,diffx,param);
    for (i=0;i<param->length;i++)
      x[i] = x[i] + h/2 * (3*diffx[i] - olddiffx[i]);
/* The following code, as a replacement for the above two lines, didn't
   seem to speed it up any.
    int l,m,c;
    for (l=0;l<=param->max_order;l+=2)
      for (m=0;m<=l;m++)
        for (c=0;c<=1;c++)
          x[ind(l,m,c)] = x[ind(l,m,c)] + h/2 * (3*diffx[ind(l,m,c)] - olddiffx[ind(l,m,c)]);
*/
    memcpy(olddiffx,diffx,param->length*sizeof(REAL));
  }
  *t += h;
}
