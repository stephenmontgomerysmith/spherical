/*
 * Created from psidot-ard.conf by expand-iterate.pl.
 */

#include "spherical.h"

extern int max_order;
static int first_mc_0 = 1;
static COMPLEX *mult_constant_0;
#define mc_count_0 39
static void initialize_mc_0();

/* Declarations of threading go here. */

extern double gamm[3][3];
extern double normgamma;
extern double b1, b2, b3, b4, b5;

void compute_psidot_ard(COMPLEX* psidot, COMPLEX* psi) {
  double a[3][3], a2[3][3], gamm2[3][3];
  int i1,i2,i3;
  double Dr[3][3];

  tensor2(psi,a);
  memset(a2,0,sizeof(a2));
  memset(gamm2,0,sizeof(gamm2));
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++) for (i3=0;i3<3;i3++) {
    a2[i1][i3] += a[i1][i2]*a[i2][i3];
    gamm2[i1][i3] += gamm[i1][i2]*gamm[i2][i3];
  }
  memset(Dr,0,sizeof(Dr));
  for (i1=0;i1<3;i1++)
    Dr[i1][i1] = b1*normgamma;
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++)
    Dr[i1][i2] += b2*normgamma*a[i1][i2] + b3*normgamma*a2[i1][i2]
                  + b4/2*gamm[i1][i2]
                  + b5/4/normgamma*gamm2[i1][i2];

  {
    if (first_mc_0) initialize_mc_0();
    /* Start-Threading */
    int lstart=0, lend=max_order+2;
    int l,m;
    COMPLEX *mc;
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      mc = mult_constant_0 + mc_count_0*ind(l,m);
      psidot[ind(l,m)] += Dr[0][0]*(mc[0]*index(psi,l-2,m-2)+mc[1]*index(psi,l-2,m)+mc[2]*index(psi,l-2,m+2)+mc[3]*index(psi,l,m-2)+mc[4]*index(psi,l,m)+mc[5]*index(psi,l,m+2)+mc[6]*index(psi,l+2,m-2)+mc[7]*index(psi,l+2,m)+mc[8]*index(psi,l+2,m+2))
                        + Dr[0][1]*(mc[9]*index(psi,l-2,m-2)+mc[10]*index(psi,l-2,m+2)+mc[11]*index(psi,l,m-2)+mc[12]*index(psi,l,m+2)+mc[13]*index(psi,l+2,m-2)+mc[14]*index(psi,l+2,m+2))
                        + Dr[0][2]*(mc[15]*index(psi,l-2,m-1)+mc[16]*index(psi,l-2,m+1)+mc[17]*index(psi,l,m-1)+mc[18]*index(psi,l,m+1)+mc[19]*index(psi,l+2,m-1)+mc[20]*index(psi,l+2,m+1))
                        + Dr[1][1]*(mc[21]*index(psi,l-2,m-2)+mc[22]*index(psi,l-2,m)+mc[23]*index(psi,l-2,m+2)+mc[24]*index(psi,l,m-2)+mc[25]*index(psi,l,m)+mc[26]*index(psi,l,m+2)+mc[27]*index(psi,l+2,m-2)+mc[28]*index(psi,l+2,m)+mc[29]*index(psi,l+2,m+2))
                        + Dr[1][2]*(mc[30]*index(psi,l-2,m-1)+mc[31]*index(psi,l-2,m+1)+mc[32]*index(psi,l,m-1)+mc[33]*index(psi,l,m+1)+mc[34]*index(psi,l+2,m-1)+mc[35]*index(psi,l+2,m+1))
                        + Dr[2][2]*(mc[36]*index(psi,l-2,m)+mc[37]*index(psi,l,m)+mc[38]*index(psi,l+2,m));
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
    if (abs(m-2)<=l-2)
      mc[0] = (((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[0] = 0;
    if (abs(m)<=l-2)
      mc[1] = (-((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[1] = 0;
    if (abs(m+2)<=l-2)
      mc[2] = (((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[2] = 0;
    if (abs(m-2)<=l)
      mc[3] = (((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[3] = 0;
    if (abs(m)<=l)
      mc[4] = ((-4*pow(ll,3)-2*pow(ll,4)-3*pow(mm,2)+pow(ll,2)*(1-2*pow(mm,2))+ll*(3-2*pow(mm,2)))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[4] = 0;
    if (abs(m+2)<=l)
      mc[5] = (((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[5] = 0;
    if (abs(m-2)<=l+2)
      mc[6] = ((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[6] = 0;
    if (abs(m)<=l+2)
      mc[7] = (-(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[7] = 0;
    if (abs(m+2)<=l+2)
      mc[8] = ((ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[8] = 0;
    if (abs(m-2)<=l-2)
      mc[9] = (0)+(-((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[9] = 0;
    if (abs(m+2)<=l-2)
      mc[10] = (0)+(((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[10] = 0;
    if (abs(m-2)<=l)
      mc[11] = (0)+(((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(6-8*ll-8*pow(ll,2)))*I;
    else
      mc[11] = 0;
    if (abs(m+2)<=l)
      mc[12] = (0)+(((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(-6+8*ll+8*pow(ll,2)))*I;
    else
      mc[12] = 0;
    if (abs(m-2)<=l+2)
      mc[13] = (0)+(-(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))*I;
    else
      mc[13] = 0;
    if (abs(m+2)<=l+2)
      mc[14] = (0)+((ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))*I;
    else
      mc[14] = 0;
    if (abs(m-1)<=l-2)
      mc[15] = (-(((-2+ll)*(1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll))))+(0)*I;
    else
      mc[15] = 0;
    if (abs(m+1)<=l-2)
      mc[16] = (((-2+ll)*(1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[16] = 0;
    if (abs(m-1)<=l)
      mc[17] = (((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[17] = 0;
    if (abs(m+1)<=l)
      mc[18] = (((3+2*ll+2*pow(ll,2))*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[18] = 0;
    if (abs(m-1)<=l+2)
      mc[19] = ((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[19] = 0;
    if (abs(m+1)<=l+2)
      mc[20] = (-((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll))))+(0)*I;
    else
      mc[20] = 0;
    if (abs(m-2)<=l-2)
      mc[21] = (-((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[21] = 0;
    if (abs(m)<=l-2)
      mc[22] = (-((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[22] = 0;
    if (abs(m+2)<=l-2)
      mc[23] = (-((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[23] = 0;
    if (abs(m-2)<=l)
      mc[24] = (-((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[24] = 0;
    if (abs(m)<=l)
      mc[25] = ((-4*pow(ll,3)-2*pow(ll,4)-3*pow(mm,2)+pow(ll,2)*(1-2*pow(mm,2))+ll*(3-2*pow(mm,2)))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[25] = 0;
    if (abs(m+2)<=l)
      mc[26] = (-((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[26] = 0;
    if (abs(m-2)<=l+2)
      mc[27] = (-(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[27] = 0;
    if (abs(m)<=l+2)
      mc[28] = (-(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[28] = 0;
    if (abs(m+2)<=l+2)
      mc[29] = (-(ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[29] = 0;
    if (abs(m-1)<=l-2)
      mc[30] = (0)+(((-2+ll)*(1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[30] = 0;
    if (abs(m+1)<=l-2)
      mc[31] = (0)+(((-2+ll)*(1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[31] = 0;
    if (abs(m-1)<=l)
      mc[32] = (0)+(((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(6-8*ll-8*pow(ll,2)))*I;
    else
      mc[32] = 0;
    if (abs(m+1)<=l)
      mc[33] = (0)+(((3+2*ll+2*pow(ll,2))*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(-6+8*ll+8*pow(ll,2)))*I;
    else
      mc[33] = 0;
    if (abs(m-1)<=l+2)
      mc[34] = (0)+(-((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll))))*I;
    else
      mc[34] = 0;
    if (abs(m+1)<=l+2)
      mc[35] = (0)+(-((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll))))*I;
    else
      mc[35] = 0;
    if (abs(m)<=l-2)
      mc[36] = (((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[36] = 0;
    if (abs(m)<=l)
      mc[37] = ((-4*pow(ll,3)-2*pow(ll,4)+3*pow(mm,2)+2*ll*pow(mm,2)+2*pow(ll,2)*(-1+pow(mm,2)))/(-3+4*ll+4*pow(ll,2)))+(0)*I;
    else
      mc[37] = 0;
    if (abs(m)<=l+2)
      mc[38] = ((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[38] = 0;
  }
}