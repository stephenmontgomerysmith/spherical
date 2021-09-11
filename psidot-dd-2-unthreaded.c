/*
 * Created from psidot-dd-2.conf by expand-iterate.pl.
 */

#include "spherical.h"

extern int max_order;
static int first_mc_0 = 1;
static COMPLEX *mult_constant_0;
#define mc_count_0 289
static void initialize_mc_0();

/* Declarations of threading go here. */

extern double gamm[3][3];
extern double normgamma;
extern double C1, C2;

void compute_psidot_dd_2(COMPLEX* psidot, COMPLEX* psi) {
  double a2[3][3], a4[3][3][3][3];

  tensor2(psi,a2);
  tensor4(psi,a4);

  {
    if (first_mc_0) initialize_mc_0();
    /* Start-Threading */
    int lstart=0, lend=max_order+2;
    int l,m;
    COMPLEX *mc;
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      mc = mult_constant_0 + mc_count_0*ind(l,m);
      COMPLEX xx, xy, xz, yy, yz, zz, xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz;
      COMPLEX temp1;
      xx = (mc[0]*index(psi,l-2,m-2)+mc[1]*index(psi,l-2,m)+mc[2]*index(psi,l-2,m+2)+mc[3]*index(psi,l,m-2)+mc[4]*index(psi,l,m)+mc[5]*index(psi,l,m+2)+mc[6]*index(psi,l+2,m-2)+mc[7]*index(psi,l+2,m)+mc[8]*index(psi,l+2,m+2));
      xy = (mc[9]*index(psi,l-2,m-2)+mc[10]*index(psi,l-2,m+2)+mc[11]*index(psi,l,m-2)+mc[12]*index(psi,l,m+2)+mc[13]*index(psi,l+2,m-2)+mc[14]*index(psi,l+2,m+2));
      xz = (mc[15]*index(psi,l-2,m-1)+mc[16]*index(psi,l-2,m+1)+mc[17]*index(psi,l,m-1)+mc[18]*index(psi,l,m+1)+mc[19]*index(psi,l+2,m-1)+mc[20]*index(psi,l+2,m+1));
      yy = (mc[21]*index(psi,l-2,m-2)+mc[22]*index(psi,l-2,m)+mc[23]*index(psi,l-2,m+2)+mc[24]*index(psi,l,m-2)+mc[25]*index(psi,l,m)+mc[26]*index(psi,l,m+2)+mc[27]*index(psi,l+2,m-2)+mc[28]*index(psi,l+2,m)+mc[29]*index(psi,l+2,m+2));
      yz = (mc[30]*index(psi,l-2,m-1)+mc[31]*index(psi,l-2,m+1)+mc[32]*index(psi,l,m-1)+mc[33]*index(psi,l,m+1)+mc[34]*index(psi,l+2,m-1)+mc[35]*index(psi,l+2,m+1));
      zz = (mc[36]*index(psi,l-2,m)+mc[37]*index(psi,l,m)+mc[38]*index(psi,l+2,m));
      xxxx = (mc[39]*index(psi,l-2,m-2)+mc[40]*index(psi,l-2,m-4)+mc[41]*index(psi,l-2,m)+mc[42]*index(psi,l-2,m+2)+mc[43]*index(psi,l-2,m+4)+mc[44]*index(psi,l-4,m-2)+mc[45]*index(psi,l-4,m-4)+mc[46]*index(psi,l-4,m)+mc[47]*index(psi,l-4,m+2)+mc[48]*index(psi,l-4,m+4)+mc[49]*index(psi,l,m-2)+mc[50]*index(psi,l,m-4)+mc[51]*index(psi,l,m)+mc[52]*index(psi,l,m+2)+mc[53]*index(psi,l,m+4)+mc[54]*index(psi,l+2,m-2)+mc[55]*index(psi,l+2,m-4)+mc[56]*index(psi,l+2,m)+mc[57]*index(psi,l+2,m+2)+mc[58]*index(psi,l+2,m+4)+mc[59]*index(psi,l+4,m-2)+mc[60]*index(psi,l+4,m-4)+mc[61]*index(psi,l+4,m)+mc[62]*index(psi,l+4,m+2)+mc[63]*index(psi,l+4,m+4));
      xxxy = (mc[64]*index(psi,l-2,m-2)+mc[65]*index(psi,l-2,m-4)+mc[66]*index(psi,l-2,m+2)+mc[67]*index(psi,l-2,m+4)+mc[68]*index(psi,l-4,m-2)+mc[69]*index(psi,l-4,m-4)+mc[70]*index(psi,l-4,m+2)+mc[71]*index(psi,l-4,m+4)+mc[72]*index(psi,l,m-2)+mc[73]*index(psi,l,m-4)+mc[74]*index(psi,l,m+2)+mc[75]*index(psi,l,m+4)+mc[76]*index(psi,l+2,m-2)+mc[77]*index(psi,l+2,m-4)+mc[78]*index(psi,l+2,m+2)+mc[79]*index(psi,l+2,m+4)+mc[80]*index(psi,l+4,m-2)+mc[81]*index(psi,l+4,m-4)+mc[82]*index(psi,l+4,m+2)+mc[83]*index(psi,l+4,m+4));
      xxxz = (mc[84]*index(psi,l-2,m-1)+mc[85]*index(psi,l-2,m-3)+mc[86]*index(psi,l-2,m+1)+mc[87]*index(psi,l-2,m+3)+mc[88]*index(psi,l-4,m-1)+mc[89]*index(psi,l-4,m-3)+mc[90]*index(psi,l-4,m+1)+mc[91]*index(psi,l-4,m+3)+mc[92]*index(psi,l,m-1)+mc[93]*index(psi,l,m-3)+mc[94]*index(psi,l,m+1)+mc[95]*index(psi,l,m+3)+mc[96]*index(psi,l+2,m-1)+mc[97]*index(psi,l+2,m-3)+mc[98]*index(psi,l+2,m+1)+mc[99]*index(psi,l+2,m+3)+mc[100]*index(psi,l+4,m-1)+mc[101]*index(psi,l+4,m-3)+mc[102]*index(psi,l+4,m+1)+mc[103]*index(psi,l+4,m+3));
      xxyy = (mc[104]*index(psi,l-2,m-4)+mc[105]*index(psi,l-2,m)+mc[106]*index(psi,l-2,m+4)+mc[107]*index(psi,l-4,m-4)+mc[108]*index(psi,l-4,m)+mc[109]*index(psi,l-4,m+4)+mc[110]*index(psi,l,m-4)+mc[111]*index(psi,l,m)+mc[112]*index(psi,l,m+4)+mc[113]*index(psi,l+2,m-4)+mc[114]*index(psi,l+2,m)+mc[115]*index(psi,l+2,m+4)+mc[116]*index(psi,l+4,m-4)+mc[117]*index(psi,l+4,m)+mc[118]*index(psi,l+4,m+4));
      xxyz = (mc[119]*index(psi,l-2,m-1)+mc[120]*index(psi,l-2,m-3)+mc[121]*index(psi,l-2,m+1)+mc[122]*index(psi,l-2,m+3)+mc[123]*index(psi,l-4,m-1)+mc[124]*index(psi,l-4,m-3)+mc[125]*index(psi,l-4,m+1)+mc[126]*index(psi,l-4,m+3)+mc[127]*index(psi,l,m-1)+mc[128]*index(psi,l,m-3)+mc[129]*index(psi,l,m+1)+mc[130]*index(psi,l,m+3)+mc[131]*index(psi,l+2,m-1)+mc[132]*index(psi,l+2,m-3)+mc[133]*index(psi,l+2,m+1)+mc[134]*index(psi,l+2,m+3)+mc[135]*index(psi,l+4,m-1)+mc[136]*index(psi,l+4,m-3)+mc[137]*index(psi,l+4,m+1)+mc[138]*index(psi,l+4,m+3));
      xxzz = (mc[139]*index(psi,l-2,m-2)+mc[140]*index(psi,l-2,m)+mc[141]*index(psi,l-2,m+2)+mc[142]*index(psi,l-4,m-2)+mc[143]*index(psi,l-4,m)+mc[144]*index(psi,l-4,m+2)+mc[145]*index(psi,l,m-2)+mc[146]*index(psi,l,m)+mc[147]*index(psi,l,m+2)+mc[148]*index(psi,l+2,m-2)+mc[149]*index(psi,l+2,m)+mc[150]*index(psi,l+2,m+2)+mc[151]*index(psi,l+4,m-2)+mc[152]*index(psi,l+4,m)+mc[153]*index(psi,l+4,m+2));
      xyyy = (mc[154]*index(psi,l-2,m-2)+mc[155]*index(psi,l-2,m-4)+mc[156]*index(psi,l-2,m+2)+mc[157]*index(psi,l-2,m+4)+mc[158]*index(psi,l-4,m-2)+mc[159]*index(psi,l-4,m-4)+mc[160]*index(psi,l-4,m+2)+mc[161]*index(psi,l-4,m+4)+mc[162]*index(psi,l,m-2)+mc[163]*index(psi,l,m-4)+mc[164]*index(psi,l,m+2)+mc[165]*index(psi,l,m+4)+mc[166]*index(psi,l+2,m-2)+mc[167]*index(psi,l+2,m-4)+mc[168]*index(psi,l+2,m+2)+mc[169]*index(psi,l+2,m+4)+mc[170]*index(psi,l+4,m-2)+mc[171]*index(psi,l+4,m-4)+mc[172]*index(psi,l+4,m+2)+mc[173]*index(psi,l+4,m+4));
      xyyz = (mc[174]*index(psi,l-2,m-1)+mc[175]*index(psi,l-2,m-3)+mc[176]*index(psi,l-2,m+1)+mc[177]*index(psi,l-2,m+3)+mc[178]*index(psi,l-4,m-1)+mc[179]*index(psi,l-4,m-3)+mc[180]*index(psi,l-4,m+1)+mc[181]*index(psi,l-4,m+3)+mc[182]*index(psi,l,m-1)+mc[183]*index(psi,l,m-3)+mc[184]*index(psi,l,m+1)+mc[185]*index(psi,l,m+3)+mc[186]*index(psi,l+2,m-1)+mc[187]*index(psi,l+2,m-3)+mc[188]*index(psi,l+2,m+1)+mc[189]*index(psi,l+2,m+3)+mc[190]*index(psi,l+4,m-1)+mc[191]*index(psi,l+4,m-3)+mc[192]*index(psi,l+4,m+1)+mc[193]*index(psi,l+4,m+3));
      xyzz = (mc[194]*index(psi,l-2,m-2)+mc[195]*index(psi,l-2,m+2)+mc[196]*index(psi,l-4,m-2)+mc[197]*index(psi,l-4,m+2)+mc[198]*index(psi,l,m-2)+mc[199]*index(psi,l,m+2)+mc[200]*index(psi,l+2,m-2)+mc[201]*index(psi,l+2,m+2)+mc[202]*index(psi,l+4,m-2)+mc[203]*index(psi,l+4,m+2));
      xzzz = (mc[204]*index(psi,l-2,m-1)+mc[205]*index(psi,l-2,m+1)+mc[206]*index(psi,l-4,m-1)+mc[207]*index(psi,l-4,m+1)+mc[208]*index(psi,l,m-1)+mc[209]*index(psi,l,m+1)+mc[210]*index(psi,l+2,m-1)+mc[211]*index(psi,l+2,m+1)+mc[212]*index(psi,l+4,m-1)+mc[213]*index(psi,l+4,m+1));
      yyyy = (mc[214]*index(psi,l-2,m-2)+mc[215]*index(psi,l-2,m-4)+mc[216]*index(psi,l-2,m)+mc[217]*index(psi,l-2,m+2)+mc[218]*index(psi,l-2,m+4)+mc[219]*index(psi,l-4,m-2)+mc[220]*index(psi,l-4,m-4)+mc[221]*index(psi,l-4,m)+mc[222]*index(psi,l-4,m+2)+mc[223]*index(psi,l-4,m+4)+mc[224]*index(psi,l,m-2)+mc[225]*index(psi,l,m-4)+mc[226]*index(psi,l,m)+mc[227]*index(psi,l,m+2)+mc[228]*index(psi,l,m+4)+mc[229]*index(psi,l+2,m-2)+mc[230]*index(psi,l+2,m-4)+mc[231]*index(psi,l+2,m)+mc[232]*index(psi,l+2,m+2)+mc[233]*index(psi,l+2,m+4)+mc[234]*index(psi,l+4,m-2)+mc[235]*index(psi,l+4,m-4)+mc[236]*index(psi,l+4,m)+mc[237]*index(psi,l+4,m+2)+mc[238]*index(psi,l+4,m+4));
      yyyz = (mc[239]*index(psi,l-2,m-1)+mc[240]*index(psi,l-2,m-3)+mc[241]*index(psi,l-2,m+1)+mc[242]*index(psi,l-2,m+3)+mc[243]*index(psi,l-4,m-1)+mc[244]*index(psi,l-4,m-3)+mc[245]*index(psi,l-4,m+1)+mc[246]*index(psi,l-4,m+3)+mc[247]*index(psi,l,m-1)+mc[248]*index(psi,l,m-3)+mc[249]*index(psi,l,m+1)+mc[250]*index(psi,l,m+3)+mc[251]*index(psi,l+2,m-1)+mc[252]*index(psi,l+2,m-3)+mc[253]*index(psi,l+2,m+1)+mc[254]*index(psi,l+2,m+3)+mc[255]*index(psi,l+4,m-1)+mc[256]*index(psi,l+4,m-3)+mc[257]*index(psi,l+4,m+1)+mc[258]*index(psi,l+4,m+3));
      yyzz = (mc[259]*index(psi,l-2,m-2)+mc[260]*index(psi,l-2,m)+mc[261]*index(psi,l-2,m+2)+mc[262]*index(psi,l-4,m-2)+mc[263]*index(psi,l-4,m)+mc[264]*index(psi,l-4,m+2)+mc[265]*index(psi,l,m-2)+mc[266]*index(psi,l,m)+mc[267]*index(psi,l,m+2)+mc[268]*index(psi,l+2,m-2)+mc[269]*index(psi,l+2,m)+mc[270]*index(psi,l+2,m+2)+mc[271]*index(psi,l+4,m-2)+mc[272]*index(psi,l+4,m)+mc[273]*index(psi,l+4,m+2));
      yzzz = (mc[274]*index(psi,l-2,m-1)+mc[275]*index(psi,l-2,m+1)+mc[276]*index(psi,l-4,m-1)+mc[277]*index(psi,l-4,m+1)+mc[278]*index(psi,l,m-1)+mc[279]*index(psi,l,m+1)+mc[280]*index(psi,l+2,m-1)+mc[281]*index(psi,l+2,m+1)+mc[282]*index(psi,l+4,m-1)+mc[283]*index(psi,l+4,m+1));
      zzzz = (mc[284]*index(psi,l-2,m)+mc[285]*index(psi,l-4,m)+mc[286]*index(psi,l,m)+mc[287]*index(psi,l+2,m)+mc[288]*index(psi,l+4,m));
      temp1 = xxzz*pow(gamm[0][1],2)*a2[0][0] - 
              2*xxyz*gamm[0][1]*gamm[0][2]*a2[0][0] + 
              xxyy*pow(gamm[0][2],2)*a2[0][0] + 
              2*xyzz*gamm[0][1]*gamm[1][1]*a2[0][0] - 
              2*xyyz*gamm[0][2]*gamm[1][1]*a2[0][0] + 
              yyzz*pow(gamm[1][1],2)*a2[0][0] - 
              2*xyyz*gamm[0][1]*gamm[1][2]*a2[0][0] + 
              2*xzzz*gamm[0][1]*gamm[1][2]*a2[0][0] + 
              2*xyyy*gamm[0][2]*gamm[1][2]*a2[0][0] - 
              2*xyzz*gamm[0][2]*gamm[1][2]*a2[0][0] - 
              2*yyyz*gamm[1][1]*gamm[1][2]*a2[0][0] + 
              2*yzzz*gamm[1][1]*gamm[1][2]*a2[0][0] + 
              yyyy*pow(gamm[1][2],2)*a2[0][0] - 
              2*yyzz*pow(gamm[1][2],2)*a2[0][0] + 
              zzzz*pow(gamm[1][2],2)*a2[0][0] - 
              2*xyzz*gamm[0][1]*gamm[2][2]*a2[0][0] + 
              2*xyyz*gamm[0][2]*gamm[2][2]*a2[0][0] - 
              2*yyzz*gamm[1][1]*gamm[2][2]*a2[0][0] + 
              2*yyyz*gamm[1][2]*gamm[2][2]*a2[0][0] - 
              2*yzzz*gamm[1][2]*gamm[2][2]*a2[0][0] + 
              yyzz*pow(gamm[2][2],2)*a2[0][0] - 
              2*xxzz*gamm[0][0]*gamm[0][1]*a2[0][1] - 
              2*xyzz*pow(gamm[0][1],2)*a2[0][1] + 
              2*xxyz*gamm[0][0]*gamm[0][2]*a2[0][1] + 
              2*xxxz*gamm[0][1]*gamm[0][2]*a2[0][1] + 
              2*xyyz*gamm[0][1]*gamm[0][2]*a2[0][1] - 
              2*xzzz*gamm[0][1]*gamm[0][2]*a2[0][1] - 
              2*xxxy*pow(gamm[0][2],2)*a2[0][1] + 
              2*xyzz*pow(gamm[0][2],2)*a2[0][1] - 
              2*xyzz*gamm[0][0]*gamm[1][1]*a2[0][1] - 
              2*yyzz*gamm[0][1]*gamm[1][1]*a2[0][1] + 
              2*xxyz*gamm[0][2]*gamm[1][1]*a2[0][1] - 
              2*yzzz*gamm[0][2]*gamm[1][1]*a2[0][1] + 
              2*xyyz*gamm[0][0]*gamm[1][2]*a2[0][1] - 
              2*xzzz*gamm[0][0]*gamm[1][2]*a2[0][1] + 
              2*xxyz*gamm[0][1]*gamm[1][2]*a2[0][1] + 
              2*yyyz*gamm[0][1]*gamm[1][2]*a2[0][1] - 
              2*yzzz*gamm[0][1]*gamm[1][2]*a2[0][1] - 
              4*xxyy*gamm[0][2]*gamm[1][2]*a2[0][1] + 
              2*xxzz*gamm[0][2]*gamm[1][2]*a2[0][1] + 
              2*yyzz*gamm[0][2]*gamm[1][2]*a2[0][1] - 
              2*zzzz*gamm[0][2]*gamm[1][2]*a2[0][1] + 
              2*xyyz*gamm[1][1]*gamm[1][2]*a2[0][1] - 
              2*xyyy*pow(gamm[1][2],2)*a2[0][1] + 
              2*xyzz*pow(gamm[1][2],2)*a2[0][1] + 
              2*xyzz*gamm[0][0]*gamm[2][2]*a2[0][1] + 
              2*xxzz*gamm[0][1]*gamm[2][2]*a2[0][1] + 
              2*yyzz*gamm[0][1]*gamm[2][2]*a2[0][1] - 
              4*xxyz*gamm[0][2]*gamm[2][2]*a2[0][1] + 
              2*yzzz*gamm[0][2]*gamm[2][2]*a2[0][1] + 
              2*xyzz*gamm[1][1]*gamm[2][2]*a2[0][1] - 
              4*xyyz*gamm[1][2]*gamm[2][2]*a2[0][1] + 
              2*xzzz*gamm[1][2]*gamm[2][2]*a2[0][1] - 
              2*xyzz*pow(gamm[2][2],2)*a2[0][1] + 
              2*xxyz*gamm[0][0]*gamm[0][1]*a2[0][2] - 
              2*xxxz*pow(gamm[0][1],2)*a2[0][2] + 
              2*xyyz*pow(gamm[0][1],2)*a2[0][2] - 
              2*xxyy*gamm[0][0]*gamm[0][2]*a2[0][2] + 
              2*xxxy*gamm[0][1]*gamm[0][2]*a2[0][2] - 
              2*xyyy*gamm[0][1]*gamm[0][2]*a2[0][2] + 
              2*xyzz*gamm[0][1]*gamm[0][2]*a2[0][2] - 
              2*xyyz*pow(gamm[0][2],2)*a2[0][2] + 
              2*xyyz*gamm[0][0]*gamm[1][1]*a2[0][2] - 
              4*xxyz*gamm[0][1]*gamm[1][1]*a2[0][2] + 
              2*yyyz*gamm[0][1]*gamm[1][1]*a2[0][2] + 
              2*xxyy*gamm[0][2]*gamm[1][1]*a2[0][2] + 
              2*yyzz*gamm[0][2]*gamm[1][1]*a2[0][2] - 
              2*xyyz*pow(gamm[1][1],2)*a2[0][2] - 
              2*xyyy*gamm[0][0]*gamm[1][2]*a2[0][2] + 
              2*xyzz*gamm[0][0]*gamm[1][2]*a2[0][2] + 
              2*xxyy*gamm[0][1]*gamm[1][2]*a2[0][2] - 
              4*xxzz*gamm[0][1]*gamm[1][2]*a2[0][2] - 
              2*yyyy*gamm[0][1]*gamm[1][2]*a2[0][2] + 
              2*yyzz*gamm[0][1]*gamm[1][2]*a2[0][2] + 
              2*xxyz*gamm[0][2]*gamm[1][2]*a2[0][2] - 
              2*yyyz*gamm[0][2]*gamm[1][2]*a2[0][2] + 
              2*yzzz*gamm[0][2]*gamm[1][2]*a2[0][2] + 
              2*xyyy*gamm[1][1]*gamm[1][2]*a2[0][2] - 
              4*xyzz*gamm[1][1]*gamm[1][2]*a2[0][2] + 
              2*xyyz*pow(gamm[1][2],2)*a2[0][2] - 
              2*xzzz*pow(gamm[1][2],2)*a2[0][2] - 
              2*xyyz*gamm[0][0]*gamm[2][2]*a2[0][2] + 
              2*xxyz*gamm[0][1]*gamm[2][2]*a2[0][2] - 
              2*yyyz*gamm[0][1]*gamm[2][2]*a2[0][2] - 
              2*yyzz*gamm[0][2]*gamm[2][2]*a2[0][2] + 
              2*xyyz*gamm[1][1]*gamm[2][2]*a2[0][2] + 
              2*xyzz*gamm[1][2]*gamm[2][2]*a2[0][2] + 
              xxzz*pow(gamm[0][0],2)*a2[1][1] + 
              2*xyzz*gamm[0][0]*gamm[0][1]*a2[1][1] + 
              yyzz*pow(gamm[0][1],2)*a2[1][1] - 
              2*xxxz*gamm[0][0]*gamm[0][2]*a2[1][1] + 
              2*xzzz*gamm[0][0]*gamm[0][2]*a2[1][1] - 
              2*xxyz*gamm[0][1]*gamm[0][2]*a2[1][1] + 
              2*yzzz*gamm[0][1]*gamm[0][2]*a2[1][1] + 
              xxxx*pow(gamm[0][2],2)*a2[1][1] - 
              2*xxzz*pow(gamm[0][2],2)*a2[1][1] + 
              zzzz*pow(gamm[0][2],2)*a2[1][1] - 
              2*xxyz*gamm[0][0]*gamm[1][2]*a2[1][1] - 
              2*xyyz*gamm[0][1]*gamm[1][2]*a2[1][1] + 
              2*xxxy*gamm[0][2]*gamm[1][2]*a2[1][1] - 
              2*xyzz*gamm[0][2]*gamm[1][2]*a2[1][1] + 
              xxyy*pow(gamm[1][2],2)*a2[1][1] - 
              2*xxzz*gamm[0][0]*gamm[2][2]*a2[1][1] - 
              2*xyzz*gamm[0][1]*gamm[2][2]*a2[1][1] + 
              2*xxxz*gamm[0][2]*gamm[2][2]*a2[1][1] - 
              2*xzzz*gamm[0][2]*gamm[2][2]*a2[1][1] + 
              2*xxyz*gamm[1][2]*gamm[2][2]*a2[1][1] + 
              xxzz*pow(gamm[2][2],2)*a2[1][1] - 
              2*xxyz*pow(gamm[0][0],2)*a2[1][2] + 
              2*xxxz*gamm[0][0]*gamm[0][1]*a2[1][2] - 
              4*xyyz*gamm[0][0]*gamm[0][1]*a2[1][2] + 
              2*xxyz*pow(gamm[0][1],2)*a2[1][2] - 
              2*yyyz*pow(gamm[0][1],2)*a2[1][2] + 
              2*xxxy*gamm[0][0]*gamm[0][2]*a2[1][2] - 
              4*xyzz*gamm[0][0]*gamm[0][2]*a2[1][2] - 
              2*xxxx*gamm[0][1]*gamm[0][2]*a2[1][2] + 
              2*xxyy*gamm[0][1]*gamm[0][2]*a2[1][2] + 
              2*xxzz*gamm[0][1]*gamm[0][2]*a2[1][2] - 
              4*yyzz*gamm[0][1]*gamm[0][2]*a2[1][2] + 
              2*xxyz*pow(gamm[0][2],2)*a2[1][2] - 
              2*yzzz*pow(gamm[0][2],2)*a2[1][2] + 
              2*xxyz*gamm[0][0]*gamm[1][1]*a2[1][2] + 
              2*xyyz*gamm[0][1]*gamm[1][1]*a2[1][2] - 
              2*xxxy*gamm[0][2]*gamm[1][1]*a2[1][2] + 
              2*xyzz*gamm[0][2]*gamm[1][1]*a2[1][2] + 
              2*xxyy*gamm[0][0]*gamm[1][2]*a2[1][2] + 
              2*xxzz*gamm[0][0]*gamm[1][2]*a2[1][2] - 
              2*xxxy*gamm[0][1]*gamm[1][2]*a2[1][2] + 
              2*xyyy*gamm[0][1]*gamm[1][2]*a2[1][2] + 
              2*xyzz*gamm[0][1]*gamm[1][2]*a2[1][2] - 
              2*xxxz*gamm[0][2]*gamm[1][2]*a2[1][2] + 
              2*xyyz*gamm[0][2]*gamm[1][2]*a2[1][2] + 
              2*xzzz*gamm[0][2]*gamm[1][2]*a2[1][2] - 
              2*xxyy*gamm[1][1]*gamm[1][2]*a2[1][2] - 
              2*xxyz*pow(gamm[1][2],2)*a2[1][2] + 
              2*xxyz*gamm[0][0]*gamm[2][2]*a2[1][2] - 
              2*xxxz*gamm[0][1]*gamm[2][2]*a2[1][2] + 
              2*xyyz*gamm[0][1]*gamm[2][2]*a2[1][2] + 
              2*xyzz*gamm[0][2]*gamm[2][2]*a2[1][2] - 
              2*xxyz*gamm[1][1]*gamm[2][2]*a2[1][2] - 
              2*xxzz*gamm[1][2]*gamm[2][2]*a2[1][2] + 
              xxyy*pow(gamm[0][0],2)*a2[2][2] - 
              2*xxxy*gamm[0][0]*gamm[0][1]*a2[2][2] + 
              2*xyyy*gamm[0][0]*gamm[0][1]*a2[2][2] + 
              xxxx*pow(gamm[0][1],2)*a2[2][2] - 
              2*xxyy*pow(gamm[0][1],2)*a2[2][2] + 
              yyyy*pow(gamm[0][1],2)*a2[2][2] + 
              2*xyyz*gamm[0][0]*gamm[0][2]*a2[2][2] - 
              2*xxyz*gamm[0][1]*gamm[0][2]*a2[2][2] + 
              2*yyyz*gamm[0][1]*gamm[0][2]*a2[2][2] + 
              yyzz*pow(gamm[0][2],2)*a2[2][2] - 
              2*xxyy*gamm[0][0]*gamm[1][1]*a2[2][2] + 
              2*xxxy*gamm[0][1]*gamm[1][1]*a2[2][2] - 
              2*xyyy*gamm[0][1]*gamm[1][1]*a2[2][2] - 
              2*xyyz*gamm[0][2]*gamm[1][1]*a2[2][2] + 
              xxyy*pow(gamm[1][1],2)*a2[2][2] - 
              2*xxyz*gamm[0][0]*gamm[1][2]*a2[2][2] + 
              2*xxxz*gamm[0][1]*gamm[1][2]*a2[2][2] - 
              2*xyyz*gamm[0][1]*gamm[1][2]*a2[2][2] - 
              2*xyzz*gamm[0][2]*gamm[1][2]*a2[2][2] + 
              2*xxyz*gamm[1][1]*gamm[1][2]*a2[2][2] + 
              xxzz*pow(gamm[1][2],2)*a2[2][2] + 
              zz*pow(gamm[0][1],2)*a4[0][0][0][0] - 
              2*yz*gamm[0][1]*gamm[0][2]*a4[0][0][0][0] + 
              yy*pow(gamm[0][2],2)*a4[0][0][0][0] - 
              2*zz*gamm[0][0]*gamm[0][1]*a4[0][0][0][1] + 
              2*yz*gamm[0][0]*gamm[0][2]*a4[0][0][0][1] + 
              2*xz*gamm[0][1]*gamm[0][2]*a4[0][0][0][1] - 
              2*xy*pow(gamm[0][2],2)*a4[0][0][0][1] + 
              2*zz*gamm[0][1]*gamm[1][1]*a4[0][0][0][1] - 
              2*yz*gamm[0][2]*gamm[1][1]*a4[0][0][0][1] - 
              2*yz*gamm[0][1]*gamm[1][2]*a4[0][0][0][1] + 
              2*yy*gamm[0][2]*gamm[1][2]*a4[0][0][0][1] + 
              2*yz*gamm[0][0]*gamm[0][1]*a4[0][0][0][2] - 
              2*xz*pow(gamm[0][1],2)*a4[0][0][0][2] - 
              2*yy*gamm[0][0]*gamm[0][2]*a4[0][0][0][2] + 
              2*xy*gamm[0][1]*gamm[0][2]*a4[0][0][0][2] + 
              2*zz*gamm[0][1]*gamm[1][2]*a4[0][0][0][2] - 
              2*yz*gamm[0][2]*gamm[1][2]*a4[0][0][0][2] - 
              2*yz*gamm[0][1]*gamm[2][2]*a4[0][0][0][2] + 
              2*yy*gamm[0][2]*gamm[2][2]*a4[0][0][0][2] + 
              zz*pow(gamm[0][0],2)*a4[0][0][1][1] - 
              2*zz*pow(gamm[0][1],2)*a4[0][0][1][1] - 
              2*xz*gamm[0][0]*gamm[0][2]*a4[0][0][1][1] + 
              2*yz*gamm[0][1]*gamm[0][2]*a4[0][0][1][1] + 
              xx*pow(gamm[0][2],2)*a4[0][0][1][1] - 
              2*zz*gamm[0][0]*gamm[1][1]*a4[0][0][1][1] + 
              2*xz*gamm[0][2]*gamm[1][1]*a4[0][0][1][1] + 
              zz*pow(gamm[1][1],2)*a4[0][0][1][1] + 
              2*yz*gamm[0][0]*gamm[1][2]*a4[0][0][1][1] + 
              2*xz*gamm[0][1]*gamm[1][2]*a4[0][0][1][1] - 
              4*xy*gamm[0][2]*gamm[1][2]*a4[0][0][1][1] - 
              2*yz*gamm[1][1]*gamm[1][2]*a4[0][0][1][1] + 
              yy*pow(gamm[1][2],2)*a4[0][0][1][1] - 
              2*yz*pow(gamm[0][0],2)*a4[0][0][1][2] + 
              2*xz*gamm[0][0]*gamm[0][1]*a4[0][0][1][2] + 
              2*yz*pow(gamm[0][1],2)*a4[0][0][1][2] + 
              2*xy*gamm[0][0]*gamm[0][2]*a4[0][0][1][2] - 
              2*xx*gamm[0][1]*gamm[0][2]*a4[0][0][1][2] - 
              2*yy*gamm[0][1]*gamm[0][2]*a4[0][0][1][2] - 
              2*zz*gamm[0][1]*gamm[0][2]*a4[0][0][1][2] + 
              2*yz*pow(gamm[0][2],2)*a4[0][0][1][2] + 
              2*yz*gamm[0][0]*gamm[1][1]*a4[0][0][1][2] - 
              4*xz*gamm[0][1]*gamm[1][1]*a4[0][0][1][2] + 
              2*xy*gamm[0][2]*gamm[1][1]*a4[0][0][1][2] - 
              2*yy*gamm[0][0]*gamm[1][2]*a4[0][0][1][2] - 
              2*zz*gamm[0][0]*gamm[1][2]*a4[0][0][1][2] + 
              2*xy*gamm[0][1]*gamm[1][2]*a4[0][0][1][2] + 
              2*xz*gamm[0][2]*gamm[1][2]*a4[0][0][1][2] + 
              2*zz*gamm[1][1]*gamm[1][2]*a4[0][0][1][2] - 
              2*yz*pow(gamm[1][2],2)*a4[0][0][1][2] + 
              2*yz*gamm[0][0]*gamm[2][2]*a4[0][0][1][2] + 
              2*xz*gamm[0][1]*gamm[2][2]*a4[0][0][1][2] - 
              4*xy*gamm[0][2]*gamm[2][2]*a4[0][0][1][2] - 
              2*yz*gamm[1][1]*gamm[2][2]*a4[0][0][1][2] + 
              2*yy*gamm[1][2]*gamm[2][2]*a4[0][0][1][2] + 
              yy*pow(gamm[0][0],2)*a4[0][0][2][2] - 
              2*xy*gamm[0][0]*gamm[0][1]*a4[0][0][2][2] + 
              xx*pow(gamm[0][1],2)*a4[0][0][2][2] + 
              2*yz*gamm[0][1]*gamm[0][2]*a4[0][0][2][2] - 
              2*yy*pow(gamm[0][2],2)*a4[0][0][2][2] + 
              2*yz*gamm[0][0]*gamm[1][2]*a4[0][0][2][2] - 
              4*xz*gamm[0][1]*gamm[1][2]*a4[0][0][2][2] + 
              2*xy*gamm[0][2]*gamm[1][2]*a4[0][0][2][2] + 
              zz*pow(gamm[1][2],2)*a4[0][0][2][2] - 
              2*yy*gamm[0][0]*gamm[2][2]*a4[0][0][2][2] + 
              2*xy*gamm[0][1]*gamm[2][2]*a4[0][0][2][2] - 
              2*yz*gamm[1][2]*gamm[2][2]*a4[0][0][2][2] + 
              yy*pow(gamm[2][2],2)*a4[0][0][2][2] + 
              2*zz*gamm[0][0]*gamm[0][1]*a4[0][1][1][1] - 
              2*xz*gamm[0][1]*gamm[0][2]*a4[0][1][1][1] - 
              2*zz*gamm[0][1]*gamm[1][1]*a4[0][1][1][1] - 
              2*xz*gamm[0][0]*gamm[1][2]*a4[0][1][1][1] + 
              2*yz*gamm[0][1]*gamm[1][2]*a4[0][1][1][1] + 
              2*xx*gamm[0][2]*gamm[1][2]*a4[0][1][1][1] + 
              2*xz*gamm[1][1]*gamm[1][2]*a4[0][1][1][1] - 
              2*xy*pow(gamm[1][2],2)*a4[0][1][1][1] - 
              4*yz*gamm[0][0]*gamm[0][1]*a4[0][1][1][2] + 
              2*xz*pow(gamm[0][1],2)*a4[0][1][1][2] + 
              2*zz*gamm[0][0]*gamm[0][2]*a4[0][1][1][2] + 
              2*xy*gamm[0][1]*gamm[0][2]*a4[0][1][1][2] - 
              2*xz*pow(gamm[0][2],2)*a4[0][1][1][2] + 
              2*xz*gamm[0][0]*gamm[1][1]*a4[0][1][1][2] + 
              2*yz*gamm[0][1]*gamm[1][1]*a4[0][1][1][2] - 
              2*xx*gamm[0][2]*gamm[1][1]*a4[0][1][1][2] - 
              2*zz*gamm[0][2]*gamm[1][1]*a4[0][1][1][2] - 
              2*xz*pow(gamm[1][1],2)*a4[0][1][1][2] + 
              2*xy*gamm[0][0]*gamm[1][2]*a4[0][1][1][2] - 
              2*xx*gamm[0][1]*gamm[1][2]*a4[0][1][1][2] - 
              2*yy*gamm[0][1]*gamm[1][2]*a4[0][1][1][2] - 
              2*zz*gamm[0][1]*gamm[1][2]*a4[0][1][1][2] + 
              2*yz*gamm[0][2]*gamm[1][2]*a4[0][1][1][2] + 
              2*xy*gamm[1][1]*gamm[1][2]*a4[0][1][1][2] + 
              2*xz*pow(gamm[1][2],2)*a4[0][1][1][2] - 
              2*xz*gamm[0][0]*gamm[2][2]*a4[0][1][1][2] + 
              2*yz*gamm[0][1]*gamm[2][2]*a4[0][1][1][2] + 
              2*xx*gamm[0][2]*gamm[2][2]*a4[0][1][1][2] + 
              2*xz*gamm[1][1]*gamm[2][2]*a4[0][1][1][2] - 
              4*xy*gamm[1][2]*gamm[2][2]*a4[0][1][1][2] + 
              2*yy*gamm[0][0]*gamm[0][1]*a4[0][1][2][2] - 
              2*xy*pow(gamm[0][1],2)*a4[0][1][2][2] - 
              4*yz*gamm[0][0]*gamm[0][2]*a4[0][1][2][2] + 
              2*xz*gamm[0][1]*gamm[0][2]*a4[0][1][2][2] + 
              2*xy*pow(gamm[0][2],2)*a4[0][1][2][2] - 
              2*xy*gamm[0][0]*gamm[1][1]*a4[0][1][2][2] + 
              2*xx*gamm[0][1]*gamm[1][1]*a4[0][1][2][2] + 
              2*yz*gamm[0][2]*gamm[1][1]*a4[0][1][2][2] + 
              2*xz*gamm[0][0]*gamm[1][2]*a4[0][1][2][2] + 
              2*yz*gamm[0][1]*gamm[1][2]*a4[0][1][2][2] - 
              2*xx*gamm[0][2]*gamm[1][2]*a4[0][1][2][2] - 
              2*yy*gamm[0][2]*gamm[1][2]*a4[0][1][2][2] - 
              2*zz*gamm[0][2]*gamm[1][2]*a4[0][1][2][2] - 
              4*xz*gamm[1][1]*gamm[1][2]*a4[0][1][2][2] + 
              2*xy*pow(gamm[1][2],2)*a4[0][1][2][2] + 
              2*xy*gamm[0][0]*gamm[2][2]*a4[0][1][2][2] - 
              2*xx*gamm[0][1]*gamm[2][2]*a4[0][1][2][2] - 
              2*yy*gamm[0][1]*gamm[2][2]*a4[0][1][2][2] + 
              2*yz*gamm[0][2]*gamm[2][2]*a4[0][1][2][2] + 
              2*xy*gamm[1][1]*gamm[2][2]*a4[0][1][2][2] + 
              2*xz*gamm[1][2]*gamm[2][2]*a4[0][1][2][2] - 
              2*xy*pow(gamm[2][2],2)*a4[0][1][2][2] + 
              2*yy*gamm[0][0]*gamm[0][2]*a4[0][2][2][2] - 
              2*xy*gamm[0][1]*gamm[0][2]*a4[0][2][2][2] - 
              2*xy*gamm[0][0]*gamm[1][2]*a4[0][2][2][2] + 
              2*xx*gamm[0][1]*gamm[1][2]*a4[0][2][2][2] + 
              2*yz*gamm[0][2]*gamm[1][2]*a4[0][2][2][2] - 
              2*xz*pow(gamm[1][2],2)*a4[0][2][2][2] - 
              2*yy*gamm[0][2]*gamm[2][2]*a4[0][2][2][2] + 
              2*xy*gamm[1][2]*gamm[2][2]*a4[0][2][2][2] + 
              zz*pow(gamm[0][1],2)*a4[1][1][1][1] - 
              2*xz*gamm[0][1]*gamm[1][2]*a4[1][1][1][1] + 
              xx*pow(gamm[1][2],2)*a4[1][1][1][1] - 
              2*yz*pow(gamm[0][1],2)*a4[1][1][1][2] + 
              2*zz*gamm[0][1]*gamm[0][2]*a4[1][1][1][2] + 
              2*xz*gamm[0][1]*gamm[1][1]*a4[1][1][1][2] + 
              2*xy*gamm[0][1]*gamm[1][2]*a4[1][1][1][2] - 
              2*xz*gamm[0][2]*gamm[1][2]*a4[1][1][1][2] - 
              2*xx*gamm[1][1]*gamm[1][2]*a4[1][1][1][2] - 
              2*xz*gamm[0][1]*gamm[2][2]*a4[1][1][1][2] + 
              2*xx*gamm[1][2]*gamm[2][2]*a4[1][1][1][2] + 
              yy*pow(gamm[0][1],2)*a4[1][1][2][2] - 
              4*yz*gamm[0][1]*gamm[0][2]*a4[1][1][2][2] + 
              zz*pow(gamm[0][2],2)*a4[1][1][2][2] - 
              2*xy*gamm[0][1]*gamm[1][1]*a4[1][1][2][2] + 
              2*xz*gamm[0][2]*gamm[1][1]*a4[1][1][2][2] + 
              xx*pow(gamm[1][1],2)*a4[1][1][2][2] + 
              2*xz*gamm[0][1]*gamm[1][2]*a4[1][1][2][2] + 
              2*xy*gamm[0][2]*gamm[1][2]*a4[1][1][2][2] - 
              2*xx*pow(gamm[1][2],2)*a4[1][1][2][2] + 
              2*xy*gamm[0][1]*gamm[2][2]*a4[1][1][2][2] - 
              2*xz*gamm[0][2]*gamm[2][2]*a4[1][1][2][2] - 
              2*xx*gamm[1][1]*gamm[2][2]*a4[1][1][2][2] + 
              xx*pow(gamm[2][2],2)*a4[1][1][2][2] + 
              2*yy*gamm[0][1]*gamm[0][2]*a4[1][2][2][2] - 
              2*yz*pow(gamm[0][2],2)*a4[1][2][2][2] - 
              2*xy*gamm[0][2]*gamm[1][1]*a4[1][2][2][2] - 
              2*xy*gamm[0][1]*gamm[1][2]*a4[1][2][2][2] + 
              2*xz*gamm[0][2]*gamm[1][2]*a4[1][2][2][2] + 
              2*xx*gamm[1][1]*gamm[1][2]*a4[1][2][2][2] + 
              2*xy*gamm[0][2]*gamm[2][2]*a4[1][2][2][2] - 
              2*xx*gamm[1][2]*gamm[2][2]*a4[1][2][2][2] + 
              yy*pow(gamm[0][2],2)*a4[2][2][2][2] - 
              2*xy*gamm[0][2]*gamm[1][2]*a4[2][2][2][2] + 
              xx*pow(gamm[1][2],2)*a4[2][2][2][2];
      psidot[ind(l,m)] += -l*(l+1)*(C1/normgamma*temp1/4 + C2*normgamma*psi[ind(l,m)]);
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
      mc[0] = (-((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll))))+(0)*I;
    else
      mc[0] = 0;
    if (abs(m)<=l-2)
      mc[1] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[1] = 0;
    if (abs(m+2)<=l-2)
      mc[2] = (-((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll))))+(0)*I;
    else
      mc[2] = 0;
    if (abs(m-2)<=l)
      mc[3] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(6-8*ll-8*pow(ll,2)))+(0)*I;
    else
      mc[3] = 0;
    if (abs(m)<=l)
      mc[4] = ((-1+ll+pow(ll,2)+pow(mm,2))/(-3+4*ll+4*pow(ll,2)))+(0)*I;
    else
      mc[4] = 0;
    if (abs(m+2)<=l)
      mc[5] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(6-8*ll-8*pow(ll,2)))+(0)*I;
    else
      mc[5] = 0;
    if (abs(m-2)<=l+2)
      mc[6] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)))+(0)*I;
    else
      mc[6] = 0;
    if (abs(m)<=l+2)
      mc[7] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll))))+(0)*I;
    else
      mc[7] = 0;
    if (abs(m+2)<=l+2)
      mc[8] = ((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)))+(0)*I;
    else
      mc[8] = 0;
    if (abs(m-2)<=l-2)
      mc[9] = (0)+((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)))*I;
    else
      mc[9] = 0;
    if (abs(m+2)<=l-2)
      mc[10] = (0)+(-((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll))))*I;
    else
      mc[10] = 0;
    if (abs(m-2)<=l)
      mc[11] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(-6+8*ll+8*pow(ll,2)))*I;
    else
      mc[11] = 0;
    if (abs(m+2)<=l)
      mc[12] = (0)+((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(6-8*ll-8*pow(ll,2)))*I;
    else
      mc[12] = 0;
    if (abs(m-2)<=l+2)
      mc[13] = (0)+(-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll))))*I;
    else
      mc[13] = 0;
    if (abs(m+2)<=l+2)
      mc[14] = (0)+((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)))*I;
    else
      mc[14] = 0;
    if (abs(m-1)<=l-2)
      mc[15] = ((sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[15] = 0;
    if (abs(m+1)<=l-2)
      mc[16] = (-((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll))))+(0)*I;
    else
      mc[16] = 0;
    if (abs(m-1)<=l)
      mc[17] = (((1-2*mm)*sqrt(1+ll-mm)*sqrt(ll+mm))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[17] = 0;
    if (abs(m+1)<=l)
      mc[18] = ((sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(6-8*ll-8*pow(ll,2)))+(0)*I;
    else
      mc[18] = 0;
    if (abs(m-1)<=l+2)
      mc[19] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)))+(0)*I;
    else
      mc[19] = 0;
    if (abs(m+1)<=l+2)
      mc[20] = (-((sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll))))+(0)*I;
    else
      mc[20] = 0;
    if (abs(m-2)<=l-2)
      mc[21] = ((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[21] = 0;
    if (abs(m)<=l-2)
      mc[22] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[22] = 0;
    if (abs(m+2)<=l-2)
      mc[23] = ((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[23] = 0;
    if (abs(m-2)<=l)
      mc[24] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[24] = 0;
    if (abs(m)<=l)
      mc[25] = ((-1+ll+pow(ll,2)+pow(mm,2))/(-3+4*ll+4*pow(ll,2)))+(0)*I;
    else
      mc[25] = 0;
    if (abs(m+2)<=l)
      mc[26] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(-6+8*ll+8*pow(ll,2)))+(0)*I;
    else
      mc[26] = 0;
    if (abs(m-2)<=l+2)
      mc[27] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll))))+(0)*I;
    else
      mc[27] = 0;
    if (abs(m)<=l+2)
      mc[28] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll))))+(0)*I;
    else
      mc[28] = 0;
    if (abs(m+2)<=l+2)
      mc[29] = (-((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll))))+(0)*I;
    else
      mc[29] = 0;
    if (abs(m-1)<=l-2)
      mc[30] = (0)+(-((sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll))))*I;
    else
      mc[30] = 0;
    if (abs(m+1)<=l-2)
      mc[31] = (0)+(-((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll))))*I;
    else
      mc[31] = 0;
    if (abs(m-1)<=l)
      mc[32] = (0)+((sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(-6+8*ll+8*pow(ll,2)))*I;
    else
      mc[32] = 0;
    if (abs(m+1)<=l)
      mc[33] = (0)+((sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(6-8*ll-8*pow(ll,2)))*I;
    else
      mc[33] = 0;
    if (abs(m-1)<=l+2)
      mc[34] = (0)+(-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll))))*I;
    else
      mc[34] = 0;
    if (abs(m+1)<=l+2)
      mc[35] = (0)+(-((sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll))))*I;
    else
      mc[35] = 0;
    if (abs(m)<=l-2)
      mc[36] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)))+(0)*I;
    else
      mc[36] = 0;
    if (abs(m)<=l)
      mc[37] = ((-1+2*ll+2*pow(ll,2)-2*pow(mm,2))/(-3+4*ll+4*pow(ll,2)))+(0)*I;
    else
      mc[37] = 0;
    if (abs(m)<=l+2)
      mc[38] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)))+(0)*I;
    else
      mc[38] = 0;
    if (abs(m-2)<=l-2)
      mc[39] = ((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-5+2*pow(ll,2)-3*mm-2*ll*mm+2*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[39] = 0;
    if (abs(m-4)<=l-2)
      mc[40] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[40] = 0;
    if (abs(m)<=l-2)
      mc[41] = ((-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-4-ll+pow(ll,2)+pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[41] = 0;
    if (abs(m+2)<=l-2)
      mc[42] = ((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(-5+2*pow(ll,2)+3*mm+2*ll*mm+2*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[42] = 0;
    if (abs(m+4)<=l-2)
      mc[43] = (-(sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[43] = 0;
    if (abs(m-2)<=l-4)
      mc[44] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3))))+(0)*I;
    else
      mc[44] = 0;
    if (abs(m-4)<=l-4)
      mc[45] = ((sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[45] = 0;
    if (abs(m)<=l-4)
      mc[46] = ((3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[46] = 0;
    if (abs(m+2)<=l-4)
      mc[47] = ((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3))))+(0)*I;
    else
      mc[47] = 0;
    if (abs(m+4)<=l-4)
      mc[48] = ((sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[48] = 0;
    if (abs(m-2)<=l)
      mc[49] = ((-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+ll+pow(ll,2)-2*mm+pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[49] = 0;
    if (abs(m-4)<=l)
      mc[50] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[50] = 0;
    if (abs(m)<=l)
      mc[51] = ((18*pow(ll,3)+9*pow(ll,4)+6*ll*(-7+pow(mm,2))+pow(ll,2)*(-33+6*pow(mm,2))+9*(4-5*pow(mm,2)+pow(mm,4)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[51] = 0;
    if (abs(m+2)<=l)
      mc[52] = ((-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+ll+pow(ll,2)+2*mm+pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[52] = 0;
    if (abs(m+4)<=l)
      mc[53] = ((3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[53] = 0;
    if (abs(m-2)<=l+2)
      mc[54] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(-3+2*pow(ll,2)-mm+2*pow(mm,2)+2*ll*(2+mm)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[54] = 0;
    if (abs(m-4)<=l+2)
      mc[55] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[55] = 0;
    if (abs(m)<=l+2)
      mc[56] = ((-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-2+3*ll+pow(ll,2)+pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[56] = 0;
    if (abs(m+2)<=l+2)
      mc[57] = ((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-3+2*pow(ll,2)-2*ll*(-2+mm)+mm+2*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[57] = 0;
    if (abs(m+4)<=l+2)
      mc[58] = (-(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[58] = 0;
    if (abs(m-2)<=l+4)
      mc[59] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)))))+(0)*I;
    else
      mc[59] = 0;
    if (abs(m-4)<=l+4)
      mc[60] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[60] = 0;
    if (abs(m)<=l+4)
      mc[61] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[61] = 0;
    if (abs(m+2)<=l+4)
      mc[62] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)))))+(0)*I;
    else
      mc[62] = 0;
    if (abs(m+4)<=l+4)
      mc[63] = ((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[63] = 0;
    if (abs(m-2)<=l-2)
      mc[64] = (0)+((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-2*pow(ll,2)+3*mm+2*ll*mm-2*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[64] = 0;
    if (abs(m-4)<=l-2)
      mc[65] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*sqrt(1+2*ll)*(60-104*ll-48*pow(ll,2)+32*pow(ll,3))))*I;
    else
      mc[65] = 0;
    if (abs(m+2)<=l-2)
      mc[66] = (0)+((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(-5+2*pow(ll,2)+3*mm+2*ll*mm+2*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[66] = 0;
    if (abs(m+4)<=l-2)
      mc[67] = (0)+(-(sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3))))*I;
    else
      mc[67] = 0;
    if (abs(m-2)<=l-4)
      mc[68] = (0)+((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[68] = 0;
    if (abs(m-4)<=l-4)
      mc[69] = (0)+((sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(240-736*ll+576*pow(ll,2)-128*pow(ll,3))))*I;
    else
      mc[69] = 0;
    if (abs(m+2)<=l-4)
      mc[70] = (0)+((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))*I;
    else
      mc[70] = 0;
    if (abs(m+4)<=l-4)
      mc[71] = (0)+((sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[71] = 0;
    if (abs(m-2)<=l)
      mc[72] = (0)+((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+ll+pow(ll,2)-2*mm+pow(mm,2)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[72] = 0;
    if (abs(m-4)<=l)
      mc[73] = (0)+((-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[73] = 0;
    if (abs(m+2)<=l)
      mc[74] = (0)+((-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+ll+pow(ll,2)+2*mm+pow(mm,2)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[74] = 0;
    if (abs(m+4)<=l)
      mc[75] = (0)+((3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[75] = 0;
    if (abs(m-2)<=l+2)
      mc[76] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(3-2*pow(ll,2)+mm-2*pow(mm,2)-2*ll*(2+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[76] = 0;
    if (abs(m-4)<=l+2)
      mc[77] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))*I;
    else
      mc[77] = 0;
    if (abs(m+2)<=l+2)
      mc[78] = (0)+((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-3+2*pow(ll,2)-2*ll*(-2+mm)+mm+2*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[78] = 0;
    if (abs(m+4)<=l+2)
      mc[79] = (0)+(-(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))*I;
    else
      mc[79] = 0;
    if (abs(m-2)<=l+4)
      mc[80] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))*I;
    else
      mc[80] = 0;
    if (abs(m-4)<=l+4)
      mc[81] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[81] = 0;
    if (abs(m+2)<=l+4)
      mc[82] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[82] = 0;
    if (abs(m+4)<=l+4)
      mc[83] = (0)+((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[83] = 0;
    if (abs(m-1)<=l-2)
      mc[84] = ((3*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(7+ll-2*pow(ll,2)+3*mm+2*ll*mm-4*pow(mm,2)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[84] = 0;
    if (abs(m-3)<=l-2)
      mc[85] = (((5+2*ll-4*mm)*sqrt(1+ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[85] = 0;
    if (abs(m+1)<=l-2)
      mc[86] = ((3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)+3*mm+4*pow(mm,2)+ll*(-1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[86] = 0;
    if (abs(m+3)<=l-2)
      mc[87] = (-(sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*(5+2*ll+4*mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[87] = 0;
    if (abs(m-1)<=l-4)
      mc[88] = ((3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[88] = 0;
    if (abs(m-3)<=l-4)
      mc[89] = ((sqrt(ll-mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))+(0)*I;
    else
      mc[89] = 0;
    if (abs(m+1)<=l-4)
      mc[90] = ((-3*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[90] = 0;
    if (abs(m+3)<=l-4)
      mc[91] = ((sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[91] = 0;
    if (abs(m-1)<=l)
      mc[92] = ((-3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-6+ll+pow(ll,2)-3*mm+3*pow(mm,2)))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[92] = 0;
    if (abs(m-3)<=l)
      mc[93] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[93] = 0;
    if (abs(m+1)<=l)
      mc[94] = ((-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(ll+pow(ll,2)+3*(-2+mm+pow(mm,2))))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[94] = 0;
    if (abs(m+3)<=l)
      mc[95] = ((3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[95] = 0;
    if (abs(m-1)<=l+2)
      mc[96] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-4+2*pow(ll,2)-mm+4*pow(mm,2)+ll*(5+2*mm)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[96] = 0;
    if (abs(m-3)<=l+2)
      mc[97] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(ll+mm)*(-3+2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[97] = 0;
    if (abs(m+1)<=l+2)
      mc[98] = ((-3*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(-4+2*pow(ll,2)+ll*(5-2*mm)+mm+4*pow(mm,2)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[98] = 0;
    if (abs(m+3)<=l+2)
      mc[99] = (((-3+2*ll-4*mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[99] = 0;
    if (abs(m-1)<=l+4)
      mc[100] = ((-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[100] = 0;
    if (abs(m-3)<=l+4)
      mc[101] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))+(0)*I;
    else
      mc[101] = 0;
    if (abs(m+1)<=l+4)
      mc[102] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[102] = 0;
    if (abs(m+3)<=l+4)
      mc[103] = (-(sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[103] = 0;
    if (abs(m-4)<=l-2)
      mc[104] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*sqrt(1+2*ll)*(60-104*ll-48*pow(ll,2)+32*pow(ll,3))))+(0)*I;
    else
      mc[104] = 0;
    if (abs(m)<=l-2)
      mc[105] = (-(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-4-ll+pow(ll,2)+pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[105] = 0;
    if (abs(m+4)<=l-2)
      mc[106] = ((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(-3+2*ll)*sqrt(1+2*ll)*(60-104*ll-48*pow(ll,2)+32*pow(ll,3))))+(0)*I;
    else
      mc[106] = 0;
    if (abs(m-4)<=l-4)
      mc[107] = ((sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(240-736*ll+576*pow(ll,2)-128*pow(ll,3))))+(0)*I;
    else
      mc[107] = 0;
    if (abs(m)<=l-4)
      mc[108] = ((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[108] = 0;
    if (abs(m+4)<=l-4)
      mc[109] = ((sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(240-736*ll+576*pow(ll,2)-128*pow(ll,3))))+(0)*I;
    else
      mc[109] = 0;
    if (abs(m-4)<=l)
      mc[110] = ((-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[110] = 0;
    if (abs(m)<=l)
      mc[111] = ((6*pow(ll,3)+3*pow(ll,4)+2*ll*(-7+pow(mm,2))+pow(ll,2)*(-11+2*pow(mm,2))+3*(4-5*pow(mm,2)+pow(mm,4)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[111] = 0;
    if (abs(m+4)<=l)
      mc[112] = ((-3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[112] = 0;
    if (abs(m-4)<=l+2)
      mc[113] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[113] = 0;
    if (abs(m)<=l+2)
      mc[114] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-2+3*ll+pow(ll,2)+pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[114] = 0;
    if (abs(m+4)<=l+2)
      mc[115] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[115] = 0;
    if (abs(m-4)<=l+4)
      mc[116] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[116] = 0;
    if (abs(m)<=l+4)
      mc[117] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))+(0)*I;
    else
      mc[117] = 0;
    if (abs(m+4)<=l+4)
      mc[118] = (-(sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[118] = 0;
    if (abs(m-1)<=l-2)
      mc[119] = (0)+((sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)-3*mm+4*pow(mm,2)-ll*(1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[119] = 0;
    if (abs(m-3)<=l-2)
      mc[120] = (0)+(-((5+2*ll-4*mm)*sqrt(1+ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[120] = 0;
    if (abs(m+1)<=l-2)
      mc[121] = (0)+((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)+3*mm+4*pow(mm,2)+ll*(-1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[121] = 0;
    if (abs(m+3)<=l-2)
      mc[122] = (0)+(-(sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*(5+2*ll+4*mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[122] = 0;
    if (abs(m-1)<=l-4)
      mc[123] = (0)+((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))*I;
    else
      mc[123] = 0;
    if (abs(m-3)<=l-4)
      mc[124] = (0)+((sqrt(ll-mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[124] = 0;
    if (abs(m+1)<=l-4)
      mc[125] = (0)+((sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))*I;
    else
      mc[125] = 0;
    if (abs(m+3)<=l-4)
      mc[126] = (0)+((sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[126] = 0;
    if (abs(m-1)<=l)
      mc[127] = (0)+((sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-6+ll+pow(ll,2)-3*mm+3*pow(mm,2)))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[127] = 0;
    if (abs(m-3)<=l)
      mc[128] = (0)+((3*(3-2*mm)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[128] = 0;
    if (abs(m+1)<=l)
      mc[129] = (0)+(-(sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(ll+pow(ll,2)+3*(-2+mm+pow(mm,2))))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[129] = 0;
    if (abs(m+3)<=l)
      mc[130] = (0)+((3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[130] = 0;
    if (abs(m-1)<=l+2)
      mc[131] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(4-2*pow(ll,2)+mm-4*pow(mm,2)-ll*(5+2*mm)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[131] = 0;
    if (abs(m-3)<=l+2)
      mc[132] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(ll+mm)*(-3+2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[132] = 0;
    if (abs(m+1)<=l+2)
      mc[133] = (0)+(-(sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(-4+2*pow(ll,2)+ll*(5-2*mm)+mm+4*pow(mm,2)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[133] = 0;
    if (abs(m+3)<=l+2)
      mc[134] = (0)+(((-3+2*ll-4*mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[134] = 0;
    if (abs(m-1)<=l+4)
      mc[135] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))*I;
    else
      mc[135] = 0;
    if (abs(m-3)<=l+4)
      mc[136] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(1+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[136] = 0;
    if (abs(m+1)<=l+4)
      mc[137] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))*I;
    else
      mc[137] = 0;
    if (abs(m+3)<=l+4)
      mc[138] = (0)+(-(sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[138] = 0;
    if (abs(m-2)<=l-2)
      mc[139] = ((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-5+4*ll*(-1+mm)+6*mm-4*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[139] = 0;
    if (abs(m)<=l-2)
      mc[140] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-1+4*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[140] = 0;
    if (abs(m+2)<=l-2)
      mc[141] = (-(sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(5+6*mm+4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[141] = 0;
    if (abs(m-2)<=l-4)
      mc[142] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[142] = 0;
    if (abs(m)<=l-4)
      mc[143] = ((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3))))+(0)*I;
    else
      mc[143] = 0;
    if (abs(m+2)<=l-4)
      mc[144] = ((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[144] = 0;
    if (abs(m-2)<=l)
      mc[145] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-6+ll+pow(ll,2)+6*mm-3*pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[145] = 0;
    if (abs(m)<=l)
      mc[146] = ((3+2*pow(ll,3)+pow(ll,4)-3*pow(mm,4)+2*pow(ll,2)*(-2+pow(mm,2))+ll*(-5+2*pow(mm,2)))/(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)))+(0)*I;
    else
      mc[146] = 0;
    if (abs(m+2)<=l)
      mc[147] = (-(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(ll+pow(ll,2)-3*(2+2*mm+pow(mm,2))))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[147] = 0;
    if (abs(m-2)<=l+2)
      mc[148] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(1+4*ll*(-1+mm)-2*mm+4*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[148] = 0;
    if (abs(m)<=l+2)
      mc[149] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-1+4*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[149] = 0;
    if (abs(m+2)<=l+2)
      mc[150] = ((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-1-2*mm-4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[150] = 0;
    if (abs(m-2)<=l+4)
      mc[151] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))))+(0)*I;
    else
      mc[151] = 0;
    if (abs(m)<=l+4)
      mc[152] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3)))))+(0)*I;
    else
      mc[152] = 0;
    if (abs(m+2)<=l+4)
      mc[153] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))))+(0)*I;
    else
      mc[153] = 0;
    if (abs(m-2)<=l-2)
      mc[154] = (0)+((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-2*pow(ll,2)+3*mm+2*ll*mm-2*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[154] = 0;
    if (abs(m-4)<=l-2)
      mc[155] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3))))*I;
    else
      mc[155] = 0;
    if (abs(m+2)<=l-2)
      mc[156] = (0)+((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(-5+2*pow(ll,2)+3*mm+2*ll*mm+2*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[156] = 0;
    if (abs(m+4)<=l-2)
      mc[157] = (0)+((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(-3+2*ll)*sqrt(1+2*ll)*(60-104*ll-48*pow(ll,2)+32*pow(ll,3))))*I;
    else
      mc[157] = 0;
    if (abs(m-2)<=l-4)
      mc[158] = (0)+((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[158] = 0;
    if (abs(m-4)<=l-4)
      mc[159] = (0)+((sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[159] = 0;
    if (abs(m+2)<=l-4)
      mc[160] = (0)+((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))*I;
    else
      mc[160] = 0;
    if (abs(m+4)<=l-4)
      mc[161] = (0)+((sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(240-736*ll+576*pow(ll,2)-128*pow(ll,3))))*I;
    else
      mc[161] = 0;
    if (abs(m-2)<=l)
      mc[162] = (0)+((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+ll+pow(ll,2)-2*mm+pow(mm,2)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[162] = 0;
    if (abs(m-4)<=l)
      mc[163] = (0)+((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[163] = 0;
    if (abs(m+2)<=l)
      mc[164] = (0)+((-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+ll+pow(ll,2)+2*mm+pow(mm,2)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[164] = 0;
    if (abs(m+4)<=l)
      mc[165] = (0)+((-3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[165] = 0;
    if (abs(m-2)<=l+2)
      mc[166] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(3-2*pow(ll,2)+mm-2*pow(mm,2)-2*ll*(2+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[166] = 0;
    if (abs(m-4)<=l+2)
      mc[167] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))*I;
    else
      mc[167] = 0;
    if (abs(m+2)<=l+2)
      mc[168] = (0)+((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-3+2*pow(ll,2)-2*ll*(-2+mm)+mm+2*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[168] = 0;
    if (abs(m+4)<=l+2)
      mc[169] = (0)+((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))*I;
    else
      mc[169] = 0;
    if (abs(m-2)<=l+4)
      mc[170] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))*I;
    else
      mc[170] = 0;
    if (abs(m-4)<=l+4)
      mc[171] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[171] = 0;
    if (abs(m+2)<=l+4)
      mc[172] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[172] = 0;
    if (abs(m+4)<=l+4)
      mc[173] = (0)+(-(sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[173] = 0;
    if (abs(m-1)<=l-2)
      mc[174] = ((sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(7+ll-2*pow(ll,2)+3*mm+2*ll*mm-4*pow(mm,2)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[174] = 0;
    if (abs(m-3)<=l-2)
      mc[175] = (-((5+2*ll-4*mm)*sqrt(1+ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[175] = 0;
    if (abs(m+1)<=l-2)
      mc[176] = ((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)+3*mm+4*pow(mm,2)+ll*(-1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[176] = 0;
    if (abs(m+3)<=l-2)
      mc[177] = ((sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*(5+2*ll+4*mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[177] = 0;
    if (abs(m-1)<=l-4)
      mc[178] = ((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[178] = 0;
    if (abs(m-3)<=l-4)
      mc[179] = ((sqrt(ll-mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[179] = 0;
    if (abs(m+1)<=l-4)
      mc[180] = ((sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))+(0)*I;
    else
      mc[180] = 0;
    if (abs(m+3)<=l-4)
      mc[181] = ((sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))+(0)*I;
    else
      mc[181] = 0;
    if (abs(m-1)<=l)
      mc[182] = (-(sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-6+ll+pow(ll,2)-3*mm+3*pow(mm,2)))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[182] = 0;
    if (abs(m-3)<=l)
      mc[183] = ((3*(3-2*mm)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[183] = 0;
    if (abs(m+1)<=l)
      mc[184] = (-(sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(ll+pow(ll,2)+3*(-2+mm+pow(mm,2))))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[184] = 0;
    if (abs(m+3)<=l)
      mc[185] = ((-3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[185] = 0;
    if (abs(m-1)<=l+2)
      mc[186] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-4+2*pow(ll,2)-mm+4*pow(mm,2)+ll*(5+2*mm)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[186] = 0;
    if (abs(m-3)<=l+2)
      mc[187] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(ll+mm)*(-3+2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[187] = 0;
    if (abs(m+1)<=l+2)
      mc[188] = (-(sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(-4+2*pow(ll,2)+ll*(5-2*mm)+mm+4*pow(mm,2)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[188] = 0;
    if (abs(m+3)<=l+2)
      mc[189] = ((sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*(3-2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[189] = 0;
    if (abs(m-1)<=l+4)
      mc[190] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[190] = 0;
    if (abs(m-3)<=l+4)
      mc[191] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(1+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[191] = 0;
    if (abs(m+1)<=l+4)
      mc[192] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))+(0)*I;
    else
      mc[192] = 0;
    if (abs(m+3)<=l+4)
      mc[193] = ((sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))+(0)*I;
    else
      mc[193] = 0;
    if (abs(m-2)<=l-2)
      mc[194] = (0)+((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-4*ll*(-1+mm)-6*mm+4*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[194] = 0;
    if (abs(m+2)<=l-2)
      mc[195] = (0)+(-(sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(5+6*mm+4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[195] = 0;
    if (abs(m-2)<=l-4)
      mc[196] = (0)+((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3))))*I;
    else
      mc[196] = 0;
    if (abs(m+2)<=l-4)
      mc[197] = (0)+((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[197] = 0;
    if (abs(m-2)<=l)
      mc[198] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-6+ll+pow(ll,2)+6*mm-3*pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[198] = 0;
    if (abs(m+2)<=l)
      mc[199] = (0)+(-(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(ll+pow(ll,2)-3*(2+2*mm+pow(mm,2))))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[199] = 0;
    if (abs(m-2)<=l+2)
      mc[200] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(1+4*ll*(-1+mm)-2*mm+4*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[200] = 0;
    if (abs(m+2)<=l+2)
      mc[201] = (0)+((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-1-2*mm-4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[201] = 0;
    if (abs(m-2)<=l+4)
      mc[202] = (0)+(-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)))))*I;
    else
      mc[202] = 0;
    if (abs(m+2)<=l+4)
      mc[203] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))))*I;
    else
      mc[203] = 0;
    if (abs(m-1)<=l-2)
      mc[204] = ((sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(8-2*pow(ll,2)+ll*(3-2*mm)-3*mm+4*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[204] = 0;
    if (abs(m+1)<=l-2)
      mc[205] = ((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-8+2*pow(ll,2)-3*mm-4*pow(mm,2)-ll*(3+2*mm)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[205] = 0;
    if (abs(m-1)<=l-4)
      mc[206] = ((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3))))+(0)*I;
    else
      mc[206] = 0;
    if (abs(m+1)<=l-4)
      mc[207] = (-((sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3)))))+(0)*I;
    else
      mc[207] = 0;
    if (abs(m-1)<=l)
      mc[208] = ((-3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-3+ll+pow(ll,2)+mm-pow(mm,2)))/(2.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[208] = 0;
    if (abs(m+1)<=l)
      mc[209] = ((-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(-3+ll+pow(ll,2)-mm-pow(mm,2)))/(2.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[209] = 0;
    if (abs(m-1)<=l+2)
      mc[210] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-3+2*pow(ll,2)+ll*(7-2*mm)+mm-4*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[210] = 0;
    if (abs(m+1)<=l+2)
      mc[211] = ((sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3-2*pow(ll,2)+mm+4*pow(mm,2)-ll*(7+2*mm)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[211] = 0;
    if (abs(m-1)<=l+4)
      mc[212] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3))))+(0)*I;
    else
      mc[212] = 0;
    if (abs(m+1)<=l+4)
      mc[213] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3)))))+(0)*I;
    else
      mc[213] = 0;
    if (abs(m-2)<=l-2)
      mc[214] = ((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-2*pow(ll,2)+3*mm+2*ll*mm-2*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[214] = 0;
    if (abs(m-4)<=l-2)
      mc[215] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[215] = 0;
    if (abs(m)<=l-2)
      mc[216] = ((-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-4-ll+pow(ll,2)+pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[216] = 0;
    if (abs(m+2)<=l-2)
      mc[217] = (-(sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(-5+2*pow(ll,2)+3*mm+2*ll*mm+2*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[217] = 0;
    if (abs(m+4)<=l-2)
      mc[218] = (-(sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[218] = 0;
    if (abs(m-2)<=l-4)
      mc[219] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[219] = 0;
    if (abs(m-4)<=l-4)
      mc[220] = ((sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[220] = 0;
    if (abs(m)<=l-4)
      mc[221] = ((3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[221] = 0;
    if (abs(m+2)<=l-4)
      mc[222] = ((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[222] = 0;
    if (abs(m+4)<=l-4)
      mc[223] = ((sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[223] = 0;
    if (abs(m-2)<=l)
      mc[224] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+ll+pow(ll,2)-2*mm+pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[224] = 0;
    if (abs(m-4)<=l)
      mc[225] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[225] = 0;
    if (abs(m)<=l)
      mc[226] = ((18*pow(ll,3)+9*pow(ll,4)+6*ll*(-7+pow(mm,2))+pow(ll,2)*(-33+6*pow(mm,2))+9*(4-5*pow(mm,2)+pow(mm,4)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[226] = 0;
    if (abs(m+2)<=l)
      mc[227] = ((3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+ll+pow(ll,2)+2*mm+pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[227] = 0;
    if (abs(m+4)<=l)
      mc[228] = ((3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[228] = 0;
    if (abs(m-2)<=l+2)
      mc[229] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(3-2*pow(ll,2)+mm-2*pow(mm,2)-2*ll*(2+mm)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[229] = 0;
    if (abs(m-4)<=l+2)
      mc[230] = (-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[230] = 0;
    if (abs(m)<=l+2)
      mc[231] = ((-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-2+3*ll+pow(ll,2)+pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[231] = 0;
    if (abs(m+2)<=l+2)
      mc[232] = (-(sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-3+2*pow(ll,2)-2*ll*(-2+mm)+mm+2*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[232] = 0;
    if (abs(m+4)<=l+2)
      mc[233] = (-(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[233] = 0;
    if (abs(m-2)<=l+4)
      mc[234] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))))+(0)*I;
    else
      mc[234] = 0;
    if (abs(m-4)<=l+4)
      mc[235] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[235] = 0;
    if (abs(m)<=l+4)
      mc[236] = ((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[236] = 0;
    if (abs(m+2)<=l+4)
      mc[237] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))))+(0)*I;
    else
      mc[237] = 0;
    if (abs(m+4)<=l+4)
      mc[238] = ((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[238] = 0;
    if (abs(m-1)<=l-2)
      mc[239] = (0)+((3*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)-3*mm+4*pow(mm,2)-ll*(1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[239] = 0;
    if (abs(m-3)<=l-2)
      mc[240] = (0)+(((5+2*ll-4*mm)*sqrt(1+ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[240] = 0;
    if (abs(m+1)<=l-2)
      mc[241] = (0)+((3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)+3*mm+4*pow(mm,2)+ll*(-1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[241] = 0;
    if (abs(m+3)<=l-2)
      mc[242] = (0)+((sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*(5+2*ll+4*mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[242] = 0;
    if (abs(m-1)<=l-4)
      mc[243] = (0)+((-3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[243] = 0;
    if (abs(m-3)<=l-4)
      mc[244] = (0)+((sqrt(ll-mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))*I;
    else
      mc[244] = 0;
    if (abs(m+1)<=l-4)
      mc[245] = (0)+((-3*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))*I;
    else
      mc[245] = 0;
    if (abs(m+3)<=l-4)
      mc[246] = (0)+((sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3))))*I;
    else
      mc[246] = 0;
    if (abs(m-1)<=l)
      mc[247] = (0)+((3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-6+ll+pow(ll,2)-3*mm+3*pow(mm,2)))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[247] = 0;
    if (abs(m-3)<=l)
      mc[248] = (0)+((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[248] = 0;
    if (abs(m+1)<=l)
      mc[249] = (0)+((-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(ll+pow(ll,2)+3*(-2+mm+pow(mm,2))))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[249] = 0;
    if (abs(m+3)<=l)
      mc[250] = (0)+((-3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))*I;
    else
      mc[250] = 0;
    if (abs(m-1)<=l+2)
      mc[251] = (0)+((-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-4+2*pow(ll,2)-mm+4*pow(mm,2)+ll*(5+2*mm)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[251] = 0;
    if (abs(m-3)<=l+2)
      mc[252] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(ll+mm)*(-3+2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[252] = 0;
    if (abs(m+1)<=l+2)
      mc[253] = (0)+((-3*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(-4+2*pow(ll,2)+ll*(5-2*mm)+mm+4*pow(mm,2)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[253] = 0;
    if (abs(m+3)<=l+2)
      mc[254] = (0)+((sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*(3-2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[254] = 0;
    if (abs(m-1)<=l+4)
      mc[255] = (0)+((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[255] = 0;
    if (abs(m-3)<=l+4)
      mc[256] = (0)+((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))*I;
    else
      mc[256] = 0;
    if (abs(m+1)<=l+4)
      mc[257] = (0)+((3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2))))*I;
    else
      mc[257] = 0;
    if (abs(m+3)<=l+4)
      mc[258] = (0)+((sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3))))*I;
    else
      mc[258] = 0;
    if (abs(m-2)<=l-2)
      mc[259] = ((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-4*ll*(-1+mm)-6*mm+4*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[259] = 0;
    if (abs(m)<=l-2)
      mc[260] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-1+4*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[260] = 0;
    if (abs(m+2)<=l-2)
      mc[261] = ((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(5+6*mm+4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[261] = 0;
    if (abs(m-2)<=l-4)
      mc[262] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3))))+(0)*I;
    else
      mc[262] = 0;
    if (abs(m)<=l-4)
      mc[263] = ((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3))))+(0)*I;
    else
      mc[263] = 0;
    if (abs(m+2)<=l-4)
      mc[264] = ((sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3))))+(0)*I;
    else
      mc[264] = 0;
    if (abs(m-2)<=l)
      mc[265] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-6+ll+pow(ll,2)+6*mm-3*pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[265] = 0;
    if (abs(m)<=l)
      mc[266] = ((3+2*pow(ll,3)+pow(ll,4)-3*pow(mm,4)+2*pow(ll,2)*(-2+pow(mm,2))+ll*(-5+2*pow(mm,2)))/(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)))+(0)*I;
    else
      mc[266] = 0;
    if (abs(m+2)<=l)
      mc[267] = ((sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(ll+pow(ll,2)-3*(2+2*mm+pow(mm,2))))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4))))+(0)*I;
    else
      mc[267] = 0;
    if (abs(m-2)<=l+2)
      mc[268] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(1+4*ll*(-1+mm)-2*mm+4*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[268] = 0;
    if (abs(m)<=l+2)
      mc[269] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-1+4*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[269] = 0;
    if (abs(m+2)<=l+2)
      mc[270] = ((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(1+2*mm+4*pow(mm,2)-4*ll*(1+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[270] = 0;
    if (abs(m-2)<=l+4)
      mc[271] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)))))+(0)*I;
    else
      mc[271] = 0;
    if (abs(m)<=l+4)
      mc[272] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3)))))+(0)*I;
    else
      mc[272] = 0;
    if (abs(m+2)<=l+4)
      mc[273] = (-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)))))+(0)*I;
    else
      mc[273] = 0;
    if (abs(m-1)<=l-2)
      mc[274] = (0)+((sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-8+2*pow(ll,2)+3*mm-4*pow(mm,2)+ll*(-3+2*mm)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[274] = 0;
    if (abs(m+1)<=l-2)
      mc[275] = (0)+((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-8+2*pow(ll,2)-3*mm-4*pow(mm,2)-ll*(3+2*mm)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))*I;
    else
      mc[275] = 0;
    if (abs(m-1)<=l-4)
      mc[276] = (0)+(-((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3)))))*I;
    else
      mc[276] = 0;
    if (abs(m+1)<=l-4)
      mc[277] = (0)+(-((sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3)))))*I;
    else
      mc[277] = 0;
    if (abs(m-1)<=l)
      mc[278] = (0)+((3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-3+ll+pow(ll,2)+mm-pow(mm,2)))/(2.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[278] = 0;
    if (abs(m+1)<=l)
      mc[279] = (0)+((-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(-3+ll+pow(ll,2)-mm-pow(mm,2)))/(2.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2))))*I;
    else
      mc[279] = 0;
    if (abs(m-1)<=l+2)
      mc[280] = (0)+(-(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-3+2*pow(ll,2)+ll*(7-2*mm)+mm-4*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[280] = 0;
    if (abs(m+1)<=l+2)
      mc[281] = (0)+((sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3-2*pow(ll,2)+mm+4*pow(mm,2)-ll*(7+2*mm)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))*I;
    else
      mc[281] = 0;
    if (abs(m-1)<=l+4)
      mc[282] = (0)+(-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3)))))*I;
    else
      mc[282] = 0;
    if (abs(m+1)<=l+4)
      mc[283] = (0)+(-((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3)))))*I;
    else
      mc[283] = 0;
    if (abs(m)<=l-2)
      mc[284] = ((2*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-7-2*ll+2*pow(ll,2)-2*pow(mm,2)))/((-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)))+(0)*I;
    else
      mc[284] = 0;
    if (abs(m)<=l-4)
      mc[285] = ((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2))))+(0)*I;
    else
      mc[285] = 0;
    if (abs(m)<=l)
      mc[286] = ((3*(3+4*pow(ll,3)+2*pow(ll,4)+10*pow(mm,2)+2*pow(mm,4)-4*ll*(2+pow(mm,2))-2*pow(ll,2)*(3+2*pow(mm,2))))/(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)))+(0)*I;
    else
      mc[286] = 0;
    if (abs(m)<=l+2)
      mc[287] = ((2*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+6*ll+2*pow(ll,2)-2*pow(mm,2)))/((-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll)))+(0)*I;
    else
      mc[287] = 0;
    if (abs(m)<=l+4)
      mc[288] = ((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(105+142*ll+60*pow(ll,2)+8*pow(ll,3))))+(0)*I;
    else
      mc[288] = 0;
  }
}
