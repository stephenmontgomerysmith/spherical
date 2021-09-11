/*
 * Created from psidot-vd.conf by expand-iterate.pl.
 */

#include "spherical.h"

extern int max_order;
static int first_mc_0 = 1;
static COMPLEX *mult_constant_0;
#define mc_count_0 1
static void initialize_mc_0();

/* Declarations of threading go here. */

extern double gamm[3][3];
extern double normgamma;
extern double C1;
extern int max_order;
extern int do_vd;

void compute_psidot_vd(COMPLEX* psidot, COMPLEX* psi) {
  int i1,i2,i3,i4;
  double a2[3][3];
  double a4[3][3][3][3];  double s;
  double Dr2;

  s = 0;
  if (do_vd == 1) {
    tensor2(psi,a2);
    for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++)
      s += a2[i1][i2]*gamm[i1][i2];
    s = fabs(s);
  } else if (do_vd ==2) {
    tensor4(psi,a4);
    for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++) for (i3=0;i3<3;i3++) for (i4=0;i4<3;i4++)
      s += a4[i1][i2][i3][i4]*gamm[i1][i2]*gamm[i3][i4];
    s = sqrt(s);
  }
  Dr2 = C1*s;

  {
    if (first_mc_0) initialize_mc_0();
    /* Start-Threading */
    int lstart=0, lend=max_order+2;
    int l,m;
    COMPLEX *mc;
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      mc = mult_constant_0 + mc_count_0*ind(l,m);
      psidot[ind(l,m)] += Dr2*(mc[0]*index(psi,l,m));
    }
    /* Stop-Threading */
  }
}

static void initialize_mc_0() {
  int l,m;
  double ll,mm;
  COMPLEX *mc;

  first_mc_0 = 0;
  mult_constant_0 = malloc(sizeof(COMPLEX)*length*mc_count_0);
  for (l=0;l<=max_order;l+=2) for (m=0;m<=l;m++) {
    ll = l;
    mm = m;
    mc = mult_constant_0 + mc_count_0*ind(l,m);
    if (abs(m)<=l)
      mc[0] = (-(ll*(1+ll)))+(0)*I;
    else
      mc[0] = 0;
  }
}
