/*
 * Created from psidot-vd.conf by expand-iterate.pl.
 */

#include "spherical.h"
static int first_mc_0 = 1;
static REAL *mult_constant_0;
#define mc_count_0 1
static void initialize_mc_0(param_list_t *param);

/* Declarations of threading go here. */

void compute_psidot_vd(REAL* psidot, REAL* psi, param_list_t *param) {
  int i1,i2;
  REAL a2[3][3];
  REAL s;
  REAL Dr2;

  s = 0;
  tensor2(psi,a2);
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++)
    s += a2[i1][i2]*param->gamm[i1*3+i2];
  Dr2 = param->normgamma*evaluate_string(param->vd_fun,'s',s/param->normgamma);

  if (first_mc_0) initialize_mc_0(param);
  {
    /* Start-Threading */
    int lstart=0, lend=param->max_order+2;
    int l,m;
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      REAL *mc = mult_constant_0 + mc_count_0*mc_ind(l,m);
      psidot[ind(l,m,0)] += Dr2*(mc[0]*psi[ind(l,m,0)]);
      psidot[ind(l,m,1)] += Dr2*(mc[0]*psi[ind(l,m,1)]);
    }
    /* Stop-Threading */
  }
}

static void initialize_mc_0(param_list_t *param) {
  int l,m;
  REAL ll,mm;
  REAL *mc;

  first_mc_0 = 0;
  mult_constant_0 = (REAL*)malloc(sizeof(REAL)*mc_ind(param->max_order+2,0)*mc_count_0);
  for (l=0;l<=param->max_order;l+=2) for (m=0;m<=l;m++) {
    ll = l;
    mm = m;
    mc = mult_constant_0 + mc_count_0*mc_ind(l,m);
    if (abs(m)<=l)
      mc[0] = -(ll*(1+ll));
    else
      mc[0] = 0;
  }
}
