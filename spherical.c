#include "spherical.h"

char Copyright[] =
"Copyright (C) 2008-2009 The Curators of the University of Missouri.\n";

int param_filename_set = 0;
const char *outfilename = "s.out";
int verbose = 0;
int verbose_print = 0;

void get_parameters(int argc, const char **argv, param_list_t *param) {
  int i;

  for (i=1;i<argc; i++) {
    if (strcmp(argv[i],"-p")==0 && i+1<argc) {
      i++;
      set_param_filename(argv[i]);
      param_filename_set = 1;
    }
    else if (strcmp(argv[i],"-o")==0 && i+1<argc) {
      i++;
      outfilename = argv[i];
    }
    else if (strcmp(argv[i],"-v")==0) {
      verbose = 1;
      if (i+1 < argc) {
        verbose_print = strtol(argv[i+1],NULL,10);
        if (verbose_print>=1) i++;
        if (verbose_print<=0) verbose_print = 0;
      }
    }
    else if (strncmp(argv[i],"-v",2)==0) {
      verbose = 1;
      verbose_print = strtol(argv[i]+2,NULL,10);
      if (verbose_print>=1) i++;
      if (verbose_print<=0) verbose_print = 0;
    }
    else {
      param_filename_set = 0;
      break;
    }
  }
  if (!param_filename_set) {
    set_param_filename("ppp");
    param_filename_set = 1;
  }
  if (!verbose) {
    verbose = 10;
    verbose_print = 10;
  }
  if (!param_filename_set) {
    printf("Usage: %s -p <parameters-file> [-o <output-file>] [-v [<number>]]\n"
           "  if -o is not set, <output-file> defaults to \"s.out\".\n"
           "  -v gives verbose output.\n"
           "  <number> is how often to print output to terminal.\n", argv[0]);
    exit(0);
  }
  set_param_verbose_level(verbose);

  param->h = param_REAL("h");
  param->print_every = param_int("print_every");
  if (getenv("NR_THREADS")==NULL)
    param->nr_threads = param_int("nr_threads");
  else {
    param->nr_threads = strtol(getenv("NR_THREADS"),NULL,10);
    if (verbose) printf("nr_threads=%d\n",param->nr_threads);
  }
  param->max_order = param_int("max_order");
/* Make data_width+2*PADDING a mulitple of 8 */
  param->data_width = (param->max_order+2*PADDING+2)/8*8;
//  param->data_width = param->max_order+2;
//  param->data_width = 56;
printf("data_width = %d\n",param->data_width);
  param->length = (param->data_width+2*8)/2*(param->data_width+2*8)*2;

  param->tstart = param_REAL("tstart");
  param->tend = param_REAL("tend");
  param->ode_adams_bash_4 = param_bool("ode_adams_bash_4");
  param->ode_rk_4 = param_bool("ode_rk_4");
  param->do_vl = param_bool("do_vl");
  if (param->do_vl) {
    param->lambda1 = param_REAL("lambda1");
    param->lambda2 = param_REAL("lambda2");
  } else
    param->lambda = param_REAL("lambda");
  param->w[0] = param_REAL("w1");
  param->w[1] = param_REAL("w2");
  param->w[2] = param_REAL("w3");
  param->gamm[0*3+0] = param_REAL("gamma11");
  param->gamm[0*3+1] = param->gamm[1*3+0] = param_REAL("gamma12");
  param->gamm[0*3+2] = param->gamm[2*3+0] = param_REAL("gamma13");
  param->gamm[1*3+1] = param_REAL("gamma22");
  param->gamm[1*3+2] = param->gamm[2*3+1] = param_REAL("gamma23");
  param->gamm[2*3+2] = param_REAL("gamma33");
  param->do_koch = param_bool("do_koch");
  param->do_dd = param_bool("do_dd");
  param->do_dd_2 = param_bool("do_dd_2");
  param->do_vd = param_bool("do_vd");
  param->do_ard = param_bool("do_ard");
  i = 0;
  if (param->do_koch) i++;
  if (param->do_dd) i++;
  if (param->do_dd_2) i++;
  if (param->do_vd) i++;
  if (param->do_ard) i++;
  if (i>1) {
    fprintf(stderr,"Error: cannot use more than one type of diffusion.\n");
    exit(1);
  }
  if (param->do_koch) {
    param->C1 = param_REAL("C1");
    param->C2 = param_REAL("C2");
    param->C3 = param_REAL("C3");
    param->CI = 0;
  } else if (param->do_dd || param->do_dd_2) {
    param->C1 = param_REAL("C1");
    param->C2 = param_REAL("C2");
    param->CI = 0;
  } else if (param->do_vd) {
    param->C1 = param_REAL("C1");
    param->CI = 0;
  } else if (param->do_ard) {
    param->b1 = param_REAL("b1");
    param->b2 = param_REAL("b2");
    param->b3 = param_REAL("b3");
    param->b4 = param_REAL("b4");
    param->b5 = param_REAL("b5");
    param->CI = 0;
  } else
    param->CI = param_REAL("CI");
  param->do_rsc = param_bool("do_rsc");
  if (param->do_rsc) param->kappa = param_REAL("kappa");
  param->print_aij = param_bool("print_aij");
  param->print_aijkl = param_bool("print_aijkl");
  param->print_aijklmn = param_bool("print_aijklmn");
  param_ignore("nr_threads");
  done_with_param();
}

int main(int argc, const char **argv) {
  int i,j;
  FILE *sout;
  REAL *psi, *psi_d;
  REAL t;
  int iteration;
  REAL a2[3][3];
  param_list_t param, *param_d;
  int do_many = 0;

  get_parameters(argc, argv, &param);
  param.normgamma = 0;
  for (i=0;i<3;i++) for (j=0;j<3;j++) param.normgamma += param.gamm[i*3+j]*param.gamm[i*3+j];
  param.normgamma = sqrt(param.normgamma/2);
  param.Dr = param.CI*param.normgamma;
  cudaMalloc((void**)&param_d,sizeof(param_list_t));
  cudaMemcpy(param_d,&param,sizeof(param_list_t),cudaMemcpyHostToDevice);

/* psi represents sqrt(4 pi) times spherical harmonic coefficients. */
  psi = (REAL*)malloc(sizeof(REAL)*param.length);
  memset(psi,0,sizeof(REAL)*param.length);
/* Start with isotropic distribution. */
  psi[ind_macro(0,0,0,param.data_width)] = 1;
  cudaMalloc((void**)&psi_d, sizeof(REAL)*param.length);
  cudaMemcpy(psi_d,psi,sizeof(REAL)*param.length,cudaMemcpyHostToDevice);

  t = param.tstart;
  iteration=0;

  sout = fopen(outfilename,"w");
  if (sout==NULL) {
    perror("unable to open output file");
    exit(1);
  }

  while (1) {
    if (iteration%param.print_every==0) {
      fprintf(sout,"%g",t);
      cudaMemcpy(psi,psi_d,sizeof(REAL)*ind_macro(4,0,0,param.data_width),cudaMemcpyDeviceToHost);
      tensor2(psi,a2,&param);
      if (!param.print_aij)
        fprintf(sout," %g %g %g %g",a2[0][0],a2[1][1],a2[2][2],a2[0][1]);
      else
        for (i=0;i<3;i++) for (j=i;j<3;j++) fprintf(sout," %g",a2[i][j]);
      if (param.print_aijkl) {
        REAL a4[3][3][3][3];
        int k,l;
        tensor4(psi,a4,&param);
        for (i=0;i<3;i++) for (j=i;j<3;j++) for (k=j;k<3;k++) for (l=k;l<3;l++)
          fprintf(sout," %g",a4[i][j][k][l]);
      }
      if (param.print_aijklmn) {
        REAL a6[3][3][3][3][3][3];
        int k,l,m,n;
        tensor6(psi,a6,&param);
        for (i=0;i<3;i++) for (j=i;j<3;j++) for (k=j;k<3;k++)
          for (l=k;l<3;l++) for (m=l;m<3;m++) for (n=m;n<3;n++)
            fprintf(sout," %g",a6[i][j][k][l][m][n]);
      }
      fprintf(sout,"\n");
      fflush(sout);
      if (verbose_print!=0 && (iteration/param.print_every)%verbose_print==0) {
        printf("%g %g %g %g %g\n",t, a2[0][0],a2[1][1],a2[2][2],a2[0][1]);
      }
    }
    if (t>=param.tend) break;
    if (iteration==param.print_every) do_many = 1;
    ode_adams_bash_2_solve(&t,psi_d,param.h,do_many,&param,param_d);
    if (do_many)
      iteration+=param.print_every;
    else
      iteration++;
  }

  exit(0);
}
