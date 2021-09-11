#include "spherical.h"

char Copyright[] =
"Copyright (C) 2008-2009 The Curators of the University of Missouri.\n";

double h;
int print_every;
int nr_threads;
int max_order;
char *outfilename = "s.out";
int verbose = 0;
int verbose_print = 0;

double tstart;
double tend;
double lambda;
double w[3];
double gamm[3][3];
int ode_adams_bash_4;
int ode_rk_4;
int do_koch;
int do_dd;
int do_dd_2;
int do_vd;
int do_ard;
double CI;
double C1;
double C2;
double C3;
double b1, b2, b3, b4, b5;
int do_vl;
double lambda1;
double lambda2;
int do_rsc;
double kappa;
int print_aij, print_aijkl, print_aijklmn;
int enter_aij_init;
double a2_init[3][3];
int param_filename_set = 0;

double Dr;
double normgamma;

void derivs(double t, COMPLEX *psi, COMPLEX *psidot) {
  t = t; /* artificial way to avoid compiler warnings */
  memset(psidot,0,sizeof(COMPLEX)*length);

  if (do_vl) compute_psidot_vl(psidot,psi);
  else compute_psidot(psidot,psi);

  if (do_koch) compute_psidot_koch(psidot,psi);
  else if (do_dd) compute_psidot_dd(psidot,psi);
  else if (do_dd_2) compute_psidot_dd_2(psidot,psi);
  else if (do_vd) compute_psidot_vd(psidot,psi);
  else if (do_ard) compute_psidot_ard(psidot,psi);

  if (do_rsc) compute_psidot_rsc(psidot,psi);
}

void get_parameters(int argc, char **argv) {
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
    printf("Usage: %s -p <parameters-file> [-o <output-file>] [-v [<number>]]\n"
           "  if -o is not set, <output-file> defaults to \"s.out\".\n"
           "  -v gives verbose output.\n"
           "  <number> is how often to print output to terminal.\n", argv[0]);
    exit(0);
  }
  set_param_verbose_level(verbose);

  h = param_double("h");
  print_every = param_int("print_every");
  if (getenv("NR_THREADS")==NULL)
    nr_threads = param_int("nr_threads");
  else {
    nr_threads = strtol(getenv("NR_THREADS"),NULL,10);
    if (verbose) printf("nr_threads=%d\n",nr_threads);
  }
  max_order = param_int("max_order");
  tstart = param_double("tstart");
  tend = param_double("tend");
  ode_adams_bash_4 = param_bool("ode_adams_bash_4");
  ode_rk_4 = param_bool("ode_rk_4");
  do_vl = param_bool("do_vl");
  if (do_vl) {
    lambda1 = param_double("lambda1");
    lambda2 = param_double("lambda2");
  } else
    lambda = param_double("lambda");
  w[0] = param_double("w1");
  w[1] = param_double("w2");
  w[2] = param_double("w3");
  gamm[0][0] = param_double("gamma11");
  gamm[0][1] = gamm[1][0] = param_double("gamma12");
  gamm[0][2] = gamm[2][0] = param_double("gamma13");
  gamm[1][1] = param_double("gamma22");
  gamm[1][2] = gamm[2][1] = param_double("gamma23");
  gamm[2][2] = param_double("gamma33");
  do_koch = param_bool("do_koch");
  do_dd = param_bool("do_dd");
  do_dd_2 = param_bool("do_dd_2");
  do_vd = param_bool("do_vd");
  do_ard = param_bool("do_ard");
  i = 0;
  if (do_koch) i++;
  if (do_dd) i++;
  if (do_dd_2) i++;
  if (do_vd) i++;
  if (do_ard) i++;
  if (i>1) {
    fprintf(stderr,"Error: cannot use more than one type of diffusion.\n");
    exit(1);
  }
  if (do_koch) {
    C1 = param_double("C1");
    C2 = param_double("C2");
    C3 = param_double("C3");
    CI = 0;
  } else if (do_dd || do_dd_2) {
    C1 = param_double("C1");
    C2 = param_double("C2");
    CI = 0;
  } else if (do_vd) {
    C1 = param_double("C1");
    CI = 0;
  } else if (do_ard) {
    b1 = param_double("b1");
    b2 = param_double("b2");
    b3 = param_double("b3");
    b4 = param_double("b4");
    b5 = param_double("b5");
    CI = 0;
  } else
    CI = param_double("CI");
  do_rsc = param_bool("do_rsc");
  if (do_rsc) kappa = param_double("kappa");
  print_aij = param_bool("print_aij");
  print_aijkl = param_bool("print_aijkl");
  print_aijklmn = param_bool("print_aijklmn");
  enter_aij_init = param_bool("enter_aij_init");
  if (enter_aij_init) {
    a2_init[0][0] = param_double("a11_init");
    a2_init[0][1] = a2_init[1][0] = param_double("a12_init");
    a2_init[0][2] = a2_init[2][0] = param_double("a13_init");
    a2_init[1][1] = param_double("a22_init");
    a2_init[1][2] = a2_init[2][1] = param_double("a23_init");
    a2_init[2][2] = param_double("a33_init");
    double eval[3], evec[3][3];
    diagonalize_sym(3,a2_init,eval,evec);
    if (fabs(eval[0]+eval[1]+eval[2] - 1.0) > 1e-15) {
      fprintf(stderr,"Warning: initial second moments tensor does not have trace equal to one.\n");
    }
    if (eval[0] < -1e-15 || eval[1] < -1e-15 || eval[2] < -1e-15) {
      fprintf(stderr,"Warning: initial second moments tensor is not positive semi-definite.\n");
    } else if (eval[0] < 0.2-1e-15 || eval[1] < 0.2-1e-15 || eval[2] < 0.2-1e-15) {
      fprintf(stderr,"Warning: initial psi is not a probability distribution because some eigenvalues are less than 1/5.\n");
    }
  }
  done_with_param();
}

int main(int argc, char **argv) {
  int i,j;
  FILE *sout;
  COMPLEX *psi;
  double t;
  int iteration;
  double a2[3][3];

  get_parameters(argc, argv);

  normgamma = 0;
  for (i=0;i<3;i++) for (j=0;j<3;j++) normgamma += gamm[i][j]*gamm[i][j];
  normgamma = sqrt(normgamma/2);
  Dr = CI*normgamma;

/* psi represents sqrt(4 pi) times spherical harmonic coefficients. */
  psi = malloc(sizeof(COMPLEX)*length);
  memset(psi,0,sizeof(COMPLEX)*length);
  if (enter_aij_init)
    reverse_tensor2(a2_init, psi);
  else
/* Start with isotropic distribution. */
    psi[0] = 1;

  t = tstart;
  iteration=0;

  sout = fopen(outfilename,"w");
  if (sout==NULL) {
    perror("unable to open output file");
    exit(1);
  }

  while (1) {
    if (iteration%print_every==0) {
      fprintf(sout,"%g",t);
      tensor2(psi,a2);
      if (!print_aij)
        fprintf(sout," %g %g %g %g",a2[0][0],a2[1][1],a2[2][2],a2[0][1]);
      else
        for (i=0;i<3;i++) for (j=i;j<3;j++) fprintf(sout," %g",a2[i][j]);
      if (print_aijkl) {
        double a4[3][3][3][3];
        int k,l;
        tensor4(psi,a4);
        for (i=0;i<3;i++) for (j=i;j<3;j++) for (k=j;k<3;k++) for (l=k;l<3;l++)
          fprintf(sout," %g",a4[i][j][k][l]);
      }
      if (print_aijklmn) {
        double a6[3][3][3][3][3][3];
        int k,l,m,n;
        tensor6(psi,a6);
        for (i=0;i<3;i++) for (j=i;j<3;j++) for (k=j;k<3;k++)
          for (l=k;l<3;l++) for (m=l;m<3;m++) for (n=m;n<3;n++)
            fprintf(sout," %g",a6[i][j][k][l][m][n]);
      }
      fprintf(sout,"\n");
      fflush(sout);
      if (verbose_print!=0 && (iteration/print_every)%verbose_print==0) {
        printf("%g %g %g %g %g\n",t, a2[0][0],a2[1][1],a2[2][2],a2[0][1]);
      }
    }
    if (t>=tend) break;
    if (ode_adams_bash_4)
      ode_adams_bash_4_solve(&t,psi,h,derivs);
    else if (ode_rk_4)
      ode_rk_4_solve(&t,psi,h,derivs);
    else
      ode_adams_bash_2_solve(&t,psi,h,derivs);
    iteration++;
  }

  exit(0);
}
