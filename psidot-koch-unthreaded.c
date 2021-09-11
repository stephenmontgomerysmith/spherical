/*
 * Created from psidot-koch.conf by expand-iterate.pl.
 */

#include "spherical.h"

extern int max_order;
static int first_mc_0 = 1;
static COMPLEX *mult_constant_0;
#define mc_count_0 53
static void initialize_mc_0();

/* Declarations of threading go here. */

extern double gamm[3][3];
extern double normgamma;
extern double C1, C2, C3;

void compute_psidot_koch(COMPLEX* psidot, COMPLEX* psi) {
  double a4[3][3][3][3], a6[3][3][3][3][3][3];
  int i1,i2,i3,i4,i5,i6;
  double D1, D2[3][3];

  tensor4(psi,a4);
  tensor6(psi,a6);
  D1 = 0;
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++)
    for (i3=0;i3<3;i3++) for (i4=0;i4<3;i4++)
      D1 += a4[i1][i2][i3][i4]*gamm[i1][i2]*gamm[i3][i4]/normgamma;
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++) {
    D2[i1][i2] = 0;
    for (i3=0;i3<3;i3++) for (i4=0;i4<3;i4++)
      for (i5=0;i5<3;i5++) for (i6=0;i6<3;i6++)
        D2[i1][i2] += a6[i1][i2][i3][i4][i5][i6]*gamm[i3][i4]*gamm[i5][i6]/normgamma;
  }

  {
    if (first_mc_0) initialize_mc_0();
    /* Start-Threading */
    int lstart=0, lend=max_order+2;
    int l,m;
    COMPLEX *mc;
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      mc = mult_constant_0 + mc_count_0*ind(l,m);
      if (C1!=0)
        psidot[ind(l,m)] += C1*D1*(mc[0]*index(psi,l,m));
      if (C2!=0)
        psidot[ind(l,m)] += C2*(D2[0][0]*(mc[1]*index(psi,l-2,m-2)+mc[2]*index(psi,l-2,m)+mc[3]*index(psi,l-2,m+2)+mc[4]*index(psi,l,m-2)+mc[5]*index(psi,l,m)+mc[6]*index(psi,l,m+2)+mc[7]*index(psi,l+2,m-2)+mc[8]*index(psi,l+2,m)+mc[9]*index(psi,l+2,m+2))
                              + D2[0][1]*(mc[10]*index(psi,l-2,m-2)+mc[11]*index(psi,l-2,m+2)+mc[12]*index(psi,l,m-2)+mc[13]*index(psi,l,m+2)+mc[14]*index(psi,l+2,m-2)+mc[15]*index(psi,l+2,m+2))
                              + D2[0][2]*(mc[16]*index(psi,l-2,m-1)+mc[17]*index(psi,l-2,m+1)+mc[18]*index(psi,l,m-1)+mc[19]*index(psi,l,m+1)+mc[20]*index(psi,l+2,m-1)+mc[21]*index(psi,l+2,m+1))
                              + D2[1][1]*(mc[22]*index(psi,l-2,m-2)+mc[23]*index(psi,l-2,m)+mc[24]*index(psi,l-2,m+2)+mc[25]*index(psi,l,m-2)+mc[26]*index(psi,l,m)+mc[27]*index(psi,l,m+2)+mc[28]*index(psi,l+2,m-2)+mc[29]*index(psi,l+2,m)+mc[30]*index(psi,l+2,m+2))
                              + D2[1][2]*(mc[31]*index(psi,l-2,m-1)+mc[32]*index(psi,l-2,m+1)+mc[33]*index(psi,l,m-1)+mc[34]*index(psi,l,m+1)+mc[35]*index(psi,l+2,m-1)+mc[36]*index(psi,l+2,m+1))
                              + D2[2][2]*(mc[37]*index(psi,l-2,m)+mc[38]*index(psi,l,m)+mc[39]*index(psi,l+2,m)));
      if (C3!=0)
        psidot[ind(l,m)] -= C3*(D2[0][0]*(mc[40]*index(psi,l,m-2)+mc[41]*index(psi,l,m)+mc[42]*index(psi,l,m+2))
                              + D2[0][1]*(mc[43]*index(psi,l,m-2)+mc[44]*index(psi,l,m+2))
                              + D2[0][2]*(mc[45]*index(psi,l,m-1)+mc[46]*index(psi,l,m+1))
                              + D2[1][1]*(mc[47]*index(psi,l,m-2)+mc[48]*index(psi,l,m)+mc[49]*index(psi,l,m+2))
                              + D2[1][2]*(mc[50]*index(psi,l,m-1)+mc[51]*index(psi,l,m+1))
                              + D2[2][2]*(mc[52]*index(psi,l,m)));
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
    if (abs(m-2)<=l-2)
      mc[1] = (((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[1] = 0;
    if (abs(m)<=l-2)
      mc[2] = (-((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[2] = 0;
    if (abs(m+2)<=l-2)
      mc[3] = (((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[3] = 0;
    if (abs(m-2)<=l)
      mc[4] = (((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[4] = 0;
    if (abs(m)<=l)
      mc[5] = ((-4*pow(ll,3)-2*pow(ll,4)-3*pow(mm,2)+pow(ll,2)*(1-2*pow(mm,2))+ll*(3-2*pow(mm,2)))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[5] = 0;
    if (abs(m+2)<=l)
      mc[6] = (((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[6] = 0;
    if (abs(m-2)<=l+2)
      mc[7] = ((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[7] = 0;
    if (abs(m)<=l+2)
      mc[8] = (-(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[8] = 0;
    if (abs(m+2)<=l+2)
      mc[9] = ((ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[9] = 0;
    if (abs(m-2)<=l-2)
      mc[10] = (0)+(-((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[10] = 0;
    if (abs(m+2)<=l-2)
      mc[11] = (0)+(((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[11] = 0;
    if (abs(m-2)<=l)
      mc[12] = (0)+(((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(6-8*ll-8*pow(ll,2)))*I;
    else
      mc[12] = 0;
    if (abs(m+2)<=l)
      mc[13] = (0)+(((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(-6+8*ll+8*pow(ll,2)))*I;
    else
      mc[13] = 0;
    if (abs(m-2)<=l+2)
      mc[14] = (0)+(-(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))*I;
    else
      mc[14] = 0;
    if (abs(m+2)<=l+2)
      mc[15] = (0)+((ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))*I;
    else
      mc[15] = 0;
    if (abs(m-1)<=l-2)
      mc[16] = (-(((-2+ll)*(1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll))))+(0)*I;
    else
      mc[16] = 0;
    if (abs(m+1)<=l-2)
      mc[17] = (((-2+ll)*(1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[17] = 0;
    if (abs(m-1)<=l)
      mc[18] = (((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[18] = 0;
    if (abs(m+1)<=l)
      mc[19] = (((3+2*ll+2*pow(ll,2))*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[19] = 0;
    if (abs(m-1)<=l+2)
      mc[20] = ((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[20] = 0;
    if (abs(m+1)<=l+2)
      mc[21] = (-((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll))))+(0)*I;
    else
      mc[21] = 0;
    if (abs(m-2)<=l-2)
      mc[22] = (-((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[22] = 0;
    if (abs(m)<=l-2)
      mc[23] = (-((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[23] = 0;
    if (abs(m+2)<=l-2)
      mc[24] = (-((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[24] = 0;
    if (abs(m-2)<=l)
      mc[25] = (-((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[25] = 0;
    if (abs(m)<=l)
      mc[26] = ((-4*pow(ll,3)-2*pow(ll,4)-3*pow(mm,2)+pow(ll,2)*(1-2*pow(mm,2))+ll*(3-2*pow(mm,2)))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[26] = 0;
    if (abs(m+2)<=l)
      mc[27] = (-((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[27] = 0;
    if (abs(m-2)<=l+2)
      mc[28] = (-(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[28] = 0;
    if (abs(m)<=l+2)
      mc[29] = (-(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[29] = 0;
    if (abs(m+2)<=l+2)
      mc[30] = (-(ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[30] = 0;
    if (abs(m-1)<=l-2)
      mc[31] = (0)+(((-2+ll)*(1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[31] = 0;
    if (abs(m+1)<=l-2)
      mc[32] = (0)+(((-2+ll)*(1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[32] = 0;
    if (abs(m-1)<=l)
      mc[33] = (0)+(((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(6-8*ll-8*pow(ll,2)))*I;
    else
      mc[33] = 0;
    if (abs(m+1)<=l)
      mc[34] = (0)+(((3+2*ll+2*pow(ll,2))*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(-6+8*ll+8*pow(ll,2)))*I;
    else
      mc[34] = 0;
    if (abs(m-1)<=l+2)
      mc[35] = (0)+(-((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll))))*I;
    else
      mc[35] = 0;
    if (abs(m+1)<=l+2)
      mc[36] = (0)+(-((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll))))*I;
    else
      mc[36] = 0;
    if (abs(m)<=l-2)
      mc[37] = (((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[37] = 0;
    if (abs(m)<=l)
      mc[38] = ((-4*pow(ll,3)-2*pow(ll,4)+3*pow(mm,2)+2*ll*pow(mm,2)+2*pow(ll,2)*(-1+pow(mm,2)))/(-3+4*ll+4*pow(ll,2)))+(0)*I;
    else
      mc[38] = 0;
    if (abs(m)<=l+2)
      mc[39] = ((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[39] = 0;
    if (abs(m-2)<=l)
      mc[40] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/4.)+(0)*I;
    else
      mc[40] = 0;
    if (abs(m)<=l)
      mc[41] = ((ll+pow(ll,2)-pow(mm,2))/2.)+(0)*I;
    else
      mc[41] = 0;
    if (abs(m+2)<=l)
      mc[42] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/4.)+(0)*I;
    else
      mc[42] = 0;
    if (abs(m-2)<=l)
      mc[43] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/2.)*I;
    else
      mc[43] = 0;
    if (abs(m+2)<=l)
      mc[44] = (0)+((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/2.)*I;
    else
      mc[44] = 0;
    if (abs(m-1)<=l)
      mc[45] = ((sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/2.)+(0)*I;
    else
      mc[45] = 0;
    if (abs(m+1)<=l)
      mc[46] = ((sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/2.)+(0)*I;
    else
      mc[46] = 0;
    if (abs(m-2)<=l)
      mc[47] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/4.)+(0)*I;
    else
      mc[47] = 0;
    if (abs(m)<=l)
      mc[48] = ((ll+pow(ll,2)-pow(mm,2))/2.)+(0)*I;
    else
      mc[48] = 0;
    if (abs(m+2)<=l)
      mc[49] = (-(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/4.)+(0)*I;
    else
      mc[49] = 0;
    if (abs(m-1)<=l)
      mc[50] = (0)+(((1-2*mm)*sqrt(1+ll-mm)*sqrt(ll+mm))/2.)*I;
    else
      mc[50] = 0;
    if (abs(m+1)<=l)
      mc[51] = (0)+((sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/2.)*I;
    else
      mc[51] = 0;
    if (abs(m)<=l)
      mc[52] = (pow(mm,2))+(0)*I;
    else
      mc[52] = 0;
  }
}
