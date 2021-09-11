#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

#define REAL double

/* PADDING must be greater than or equal to degree of F^* (see the paper).
   Also, it must be a multiple of 2. */
#ifndef PADDING
#define PADDING 4
#endif

/*
 * psi[ind(l,m,c)] for (c=0,1) represent the real and imaginary components of
 * the Y_l^m spherical harmonic coefficient of psi.
 */
#define ind(l,m,c) ((((l)+5*PADDING)*((l)+PADDING)/4+(m)+PADDING)*2+c)
#define mc_ind(l,m) ((l)*(l)/4+(m))

typedef struct {
  int verbose;
  int verbose_print;
  const char *outfilename;
  REAL h;
  int print_every;
  int nr_threads;
  int max_order;
  REAL tstart;
  REAL tend;
  REAL lambda;
  REAL w[3];
  REAL gamm[9];
  int ode_adams_bash_2;
  int ode_adams_bash_4;
  int ode_rk_4;
  int ode_rkf_45;
  REAL tol;
  int do_koch;
  int do_dd;
  int do_vd;
  char *vd_fun;
  int do_ard;
  REAL CI;
  REAL C1;
  REAL C2;
  REAL C3;
  REAL b1, b2, b3, b4, b5;
  int do_vl;
  char *vl_fun;
  int do_rsc;
  REAL kappa;
  int print_aij, print_aijkl, print_aijklmn;
  int print_daij, print_daijkl, print_daijklmn;
  REAL normgamma;
  int length;
} param_list_t;

void derivs(REAL t, REAL *psi, REAL *psidot,param_list_t *param);
void condon_shortley(REAL *psi, param_list_t *param);
void ode_adams_bash_2_solve(REAL *t, REAL *x, REAL h,
                            void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                            param_list_t *param);
void ode_adams_bash_4_solve(REAL *t, REAL *x, REAL h,
                            void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                            param_list_t *param);
void ode_rk_4_solve(REAL *t, REAL *x, REAL h,
                    void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                    param_list_t *param);
void ode_rkf_23_solve(REAL *t, REAL *x, REAL *h_use,
                      void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                      param_list_t *param);
void ode_rkf_45_solve(REAL *t, REAL *x, REAL *h_use,
                      void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param),
                      param_list_t *param);
void compute_psidot(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_koch(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_dd(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_vd(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_vl(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_ard(REAL* psidot, REAL* psi, param_list_t *param);
void tensor2(REAL *psi, REAL a[3][3]);
void tensor4(REAL *psi, REAL a[3][3][3][3]);
void tensor6(REAL *psi, REAL a[3][3][3][3][3][3]);
void compute_psidot_rsc(REAL* psidot, REAL* psi, param_list_t *param);
void reverse_tensor2(REAL a[3][3], REAL *psi);
void reverse_tensor4(REAL a[3][3][3][3], REAL *psi);
void diagonalize_sym(int n, REAL *A, REAL *eval, REAL *evec);
void get_parameters(int argc, const char **argv, param_list_t *param);
int param_bool(const char *p);
int param_int(const char *p);
REAL param_REAL(const char *p);
int param_choice(const char *p, ...);
char * param_string(const char *p);
REAL evaluate_string(char *s, char var, REAL x);
void set_param_filename(const char *f);
void set_param_verbose_level(int v);
void done_with_param();
void param_ignore(const char *p);
int check_param(const char *p);
