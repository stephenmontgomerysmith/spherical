#include "spherical.h"

void derivs(REAL t, REAL *psi, REAL *psidot,param_list_t *param) {

  memset(psidot,0,sizeof(REAL)*param->length);

  condon_shortley(psi,param);

  if (param->do_vl) compute_psidot_vl(psidot,psi,param);
  else compute_psidot(psidot,psi,param);

  if (param->do_koch) compute_psidot_koch(psidot,psi,param);
  else if (param->do_dd) compute_psidot_dd(psidot,psi,param);
  else if (param->do_vd) compute_psidot_vd(psidot,psi,param);
  else if (param->do_ard) compute_psidot_ard(psidot,psi,param);

  if (param->do_rsc) compute_psidot_rsc(psidot,psi,param);
}

/* Y_l^(-m) = (-1)^m conj(Y_l^m) */

void condon_shortley(REAL *psi, param_list_t *param) {
  int l,m,c;

  for (l=0;l<=param->max_order;l+=2)
    for (m=-PADDING;m<0;m++)
      for (c=0;c<=1;c++)
        psi[ind(l,m,c)] = (m&1)^c ? -psi[ind(l,-m,c)] : psi[ind(l,-m,c)];
}
