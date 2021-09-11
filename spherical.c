#include "spherical.h"

char Copyright[] =
"Copyright (C) 2008-2009 The Curators of the University of Missouri.\n";

int main(int argc, const char **argv) {
  int i,j;
  FILE *sout;
  REAL *psi;
  REAL t;
  int iteration;
  REAL a2[3][3];
  param_list_t param;
  REAL *psidot;

  get_parameters(argc, argv, &param);

  param.normgamma = 0;
  for (i=0;i<3;i++) for (j=0;j<3;j++) param.normgamma += param.gamm[i*3+j]*param.gamm[i*3+j];
  param.normgamma = sqrt(param.normgamma/2);

/* psi represents sqrt(4 pi) times spherical harmonic coefficients. */
  psi = (REAL*)malloc(sizeof(REAL)*param.length);
  memset(psi,0,sizeof(REAL)*param.length);
/* Start with isotropic distribution. */
  psi[ind(0,0,0)] = 1;

  if (param.print_daij || param.print_daijkl || param.print_daijklmn)
    psidot = (REAL*)malloc(sizeof(REAL)*param.length);

  t = param.tstart;
  iteration=0;

  sout = fopen(param.outfilename,"w");
  if (sout==NULL) {
    perror("unable to open output file");
    exit(1);
  }

  while (1) {
    if (iteration%param.print_every==0) {
      fprintf(sout,"%g",t);
      tensor2(psi,a2);
      if (!param.print_aij)
        fprintf(sout," %g %g %g %g",a2[0][0],a2[1][1],a2[2][2],a2[0][1]);
      else
        for (i=0;i<3;i++) for (j=i;j<3;j++) fprintf(sout," %g",a2[i][j]);
      if (param.print_aijkl) {
        REAL a4[3][3][3][3];
        int k,l;
        tensor4(psi,a4);
        for (i=0;i<3;i++) for (j=i;j<3;j++) for (k=j;k<3;k++) for (l=k;l<3;l++)
          fprintf(sout," %g",a4[i][j][k][l]);
      }
      if (param.print_aijklmn) {
        REAL a6[3][3][3][3][3][3];
        int k,l,m,n;
        tensor6(psi,a6);
        for (i=0;i<3;i++) for (j=i;j<3;j++) for (k=j;k<3;k++)
          for (l=k;l<3;l++) for (m=l;m<3;m++) for (n=m;n<3;n++)
            fprintf(sout," %g",a6[i][j][k][l][m][n]);
      }
      if (param.print_daij || param.print_daijkl || param.print_daijklmn) {
        REAL da2[3][3];
        REAL da4[3][3][3][3];
        REAL da6[3][3][3][3][3][3];
        int k,l,m,n;
        derivs(t, psi, psidot, &param);
        if (param.print_daij) {
          tensor2(psidot,da2);
          for (i=0;i<3;i++) for (j=i;j<3;j++) fprintf(sout," %g",da2[i][j]);
        }
        if (param.print_daijkl) {
          tensor4(psidot,da4);
          for (i=0;i<3;i++) for (j=i;j<3;j++) for (k=j;k<3;k++) for (l=k;l<3;l++)
            fprintf(sout," %g",da4[i][j][k][l]);
        }
        if (param.print_daijklmn) {
          tensor6(psidot,da6);
          for (i=0;i<3;i++) for (j=i;j<3;j++) for (k=j;k<3;k++)
            for (l=k;l<3;l++) for (m=l;m<3;m++) for (n=m;n<3;n++)
              fprintf(sout," %g",da6[i][j][k][l][m][n]);
        }
      }
      fprintf(sout,"\n");
      fflush(sout);
      if (param.verbose_print!=0 && (iteration/param.print_every)%param.verbose_print==0) {
        printf("%g %g %g %g %g\n",t, a2[0][0],a2[1][1],a2[2][2],a2[0][1]);
      }
    }
    if (t>=param.tend) break;
    if (param.ode_adams_bash_2)
      ode_adams_bash_2_solve(&t,psi,param.h,derivs,&param);
    else if (param.ode_adams_bash_4)
      ode_adams_bash_4_solve(&t,psi,param.h,derivs,&param);
    else if (param.ode_rk_4)
      ode_rk_4_solve(&t,psi,param.h,derivs,&param);
    else if (param.ode_rkf_45)
      ode_rkf_45_solve(&t,psi,&(param.h),derivs,&param);
    else
      ode_rkf_23_solve(&t,psi,&(param.h),derivs,&param);

    if (t+param.h*(1+1e-7)>param.tend) param.h = (1+1e-7)*(param.tend-t);
    iteration++;
  }

  if (param.verbose) printf("Total number of iterations = %d\n",iteration);
  exit(0);
}
