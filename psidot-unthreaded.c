/*
 * Created from psidot.conf by expand-iterate.pl.
 */

#include "spherical.h"

extern int max_order;
static int first_mc_0 = 1;
static COMPLEX *mult_constant_0;
#define mc_count_0 45
static void initialize_mc_0();

/* Declarations of threading go here. */

extern double w[3];
extern double gamm[3][3];
extern double lambda;
extern double Dr;

void compute_psidot(COMPLEX* psidot, COMPLEX* psi) {
  {
    if (first_mc_0) initialize_mc_0();
    /* Start-Threading */
    int lstart=0, lend=max_order+2;
    int l,m;
    COMPLEX *mc;
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      mc = mult_constant_0 + mc_count_0*ind(l,m);
  /* Jeffery's Equation */
      psidot[ind(l,m)] += w[0]*(mc[0]*index(psi,l,m-1)+mc[1]*index(psi,l,m+1))
                        + w[1]*(mc[2]*index(psi,l,m-1)+mc[3]*index(psi,l,m+1))
                        + w[2]*(mc[4]*index(psi,l,m))
                        + lambda*gamm[0][0]*(mc[5]*index(psi,l-2,m-2)+mc[6]*index(psi,l-2,m)+mc[7]*index(psi,l-2,m+2)+mc[8]*index(psi,l,m-2)+mc[9]*index(psi,l,m)+mc[10]*index(psi,l,m+2)+mc[11]*index(psi,l+2,m-2)+mc[12]*index(psi,l+2,m)+mc[13]*index(psi,l+2,m+2))
                        + lambda*gamm[0][1]*(mc[14]*index(psi,l-2,m-2)+mc[15]*index(psi,l-2,m+2)+mc[16]*index(psi,l,m-2)+mc[17]*index(psi,l,m+2)+mc[18]*index(psi,l+2,m-2)+mc[19]*index(psi,l+2,m+2))
                        + lambda*gamm[0][2]*(mc[20]*index(psi,l-2,m-1)+mc[21]*index(psi,l-2,m+1)+mc[22]*index(psi,l,m-1)+mc[23]*index(psi,l,m+1)+mc[24]*index(psi,l+2,m-1)+mc[25]*index(psi,l+2,m+1))
                        + lambda*gamm[1][1]*(mc[26]*index(psi,l-2,m-2)+mc[27]*index(psi,l-2,m)+mc[28]*index(psi,l-2,m+2)+mc[29]*index(psi,l,m-2)+mc[30]*index(psi,l,m)+mc[31]*index(psi,l,m+2)+mc[32]*index(psi,l+2,m-2)+mc[33]*index(psi,l+2,m)+mc[34]*index(psi,l+2,m+2))
                        + lambda*gamm[1][2]*(mc[35]*index(psi,l-2,m-1)+mc[36]*index(psi,l-2,m+1)+mc[37]*index(psi,l,m-1)+mc[38]*index(psi,l,m+1)+mc[39]*index(psi,l+2,m-1)+mc[40]*index(psi,l+2,m+1))
                        + lambda*gamm[2][2]*(mc[41]*index(psi,l-2,m)+mc[42]*index(psi,l,m)+mc[43]*index(psi,l+2,m))
  /* Folgar-Tucker correction */
                        + Dr*(mc[44]*index(psi,l,m));
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
    if (abs(m-1)<=l)
      mc[0] = (0)+(-(sqrt(1+ll-mm)*sqrt(ll+mm))/4.)*I;
    else
      mc[0] = 0;
    if (abs(m+1)<=l)
      mc[1] = (0)+(-(sqrt(ll-mm)*sqrt(1+ll+mm))/4.)*I;
    else
      mc[1] = 0;
    if (abs(m-1)<=l)
      mc[2] = (-(sqrt(1+ll-mm)*sqrt(ll+mm))/4.)+(0)*I;
    else
      mc[2] = 0;
    if (abs(m+1)<=l)
      mc[3] = ((sqrt(ll-mm)*sqrt(1+ll+mm))/4.)+(0)*I;
    else
      mc[3] = 0;
    if (abs(m)<=l)
      mc[4] = (0)+(-mm/2.)*I;
    else
      mc[4] = 0;
    if (abs(m-2)<=l-2)
      mc[5] = (((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[5] = 0;
    if (abs(m)<=l-2)
      mc[6] = (-((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[6] = 0;
    if (abs(m+2)<=l-2)
      mc[7] = (((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[7] = 0;
    if (abs(m-2)<=l)
      mc[8] = ((-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[8] = 0;
    if (abs(m)<=l)
      mc[9] = ((ll+pow(ll,2)-3*pow(mm,2))/(12-16*ll-16*pow(ll,2)))+(0)*I;
    else
      mc[9] = 0;
    if (abs(m+2)<=l)
      mc[10] = ((-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(8.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[10] = 0;
    if (abs(m-2)<=l+2)
      mc[11] = (-((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll))))+(0)*I;
    else
      mc[11] = 0;
    if (abs(m)<=l+2)
      mc[12] = ((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)))+(0)*I;
    else
      mc[12] = 0;
    if (abs(m+2)<=l+2)
      mc[13] = (-((ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll))))+(0)*I;
    else
      mc[13] = 0;
    if (abs(m-2)<=l-2)
      mc[14] = (0)+(-((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[14] = 0;
    if (abs(m+2)<=l-2)
      mc[15] = (0)+(((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[15] = 0;
    if (abs(m-2)<=l)
      mc[16] = (0)+((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[16] = 0;
    if (abs(m+2)<=l)
      mc[17] = (0)+((-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[17] = 0;
    if (abs(m-2)<=l+2)
      mc[18] = (0)+((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)))*I;
    else
      mc[18] = 0;
    if (abs(m+2)<=l+2)
      mc[19] = (0)+(-((ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll))))*I;
    else
      mc[19] = 0;
    if (abs(m-1)<=l-2)
      mc[20] = (-((1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[20] = 0;
    if (abs(m+1)<=l-2)
      mc[21] = (((1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[21] = 0;
    if (abs(m-1)<=l)
      mc[22] = ((3*(1-2*mm)*sqrt(1+ll-mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[22] = 0;
    if (abs(m+1)<=l)
      mc[23] = ((-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[23] = 0;
    if (abs(m-1)<=l+2)
      mc[24] = (-((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll))))+(0)*I;
    else
      mc[24] = 0;
    if (abs(m+1)<=l+2)
      mc[25] = ((ll*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)))+(0)*I;
    else
      mc[25] = 0;
    if (abs(m-2)<=l-2)
      mc[26] = (-((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[26] = 0;
    if (abs(m)<=l-2)
      mc[27] = (-((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[27] = 0;
    if (abs(m+2)<=l-2)
      mc[28] = (-((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[28] = 0;
    if (abs(m-2)<=l)
      mc[29] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[29] = 0;
    if (abs(m)<=l)
      mc[30] = ((ll+pow(ll,2)-3*pow(mm,2))/(12-16*ll-16*pow(ll,2)))+(0)*I;
    else
      mc[30] = 0;
    if (abs(m+2)<=l)
      mc[31] = ((3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(8.*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[31] = 0;
    if (abs(m-2)<=l+2)
      mc[32] = ((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll)))+(0)*I;
    else
      mc[32] = 0;
    if (abs(m)<=l+2)
      mc[33] = ((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)))+(0)*I;
    else
      mc[33] = 0;
    if (abs(m+2)<=l+2)
      mc[34] = ((ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll)))+(0)*I;
    else
      mc[34] = 0;
    if (abs(m-1)<=l-2)
      mc[35] = (0)+(((1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[35] = 0;
    if (abs(m+1)<=l-2)
      mc[36] = (0)+(((1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[36] = 0;
    if (abs(m-1)<=l)
      mc[37] = (0)+((3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[37] = 0;
    if (abs(m+1)<=l)
      mc[38] = (0)+((-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[38] = 0;
    if (abs(m-1)<=l+2)
      mc[39] = (0)+((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)))*I;
    else
      mc[39] = 0;
    if (abs(m+1)<=l+2)
      mc[40] = (0)+((ll*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)))*I;
    else
      mc[40] = 0;
    if (abs(m)<=l-2)
      mc[41] = (((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[41] = 0;
    if (abs(m)<=l)
      mc[42] = ((ll+pow(ll,2)-3*pow(mm,2))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[42] = 0;
    if (abs(m)<=l+2)
      mc[43] = (-((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll))))+(0)*I;
    else
      mc[43] = 0;
    if (abs(m)<=l)
      mc[44] = (-(ll*(1+ll)))+(0)*I;
    else
      mc[44] = 0;
  }
}
