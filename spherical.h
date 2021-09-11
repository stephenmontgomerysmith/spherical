#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef NO_STDIO
#include <stdio.h>
#endif
#include <stdarg.h>
#include <ctype.h>
#include <cuda_runtime_api.h>

#define REAL float

#ifndef PADDING
#define PADDING 4
#endif

/*
 * psi[ind(l,m,c)] for (c=0,1) represent the real and imaginary components of
 * the Y_l^m spherical harmonic coefficient of psi.
 */
#define ind_macro(l,m,c,w) \
  ((((l)+8)/2*(w+2*8)+(m)+8)*2+c)
#define ind(l,m,c) ind_macro(l,m,c,param->data_width)

typedef struct {
  REAL h;
  int print_every;
  int nr_threads;
  int max_order;
  int data_width;
  REAL tstart;
  REAL tend;
  REAL lambda;
  REAL w[3];
  REAL gamm[9];
  int ode_adams_bash_4;
  int ode_rk_4;
  int do_koch;
  int do_dd;
  int do_dd_2;
  int do_vd;
  int do_ard;
  REAL CI;
  REAL C1;
  REAL C2;
  REAL C3;
  REAL b1, b2, b3, b4, b5;
  int do_vl;
  REAL lambda1;
  REAL lambda2;
  int do_rsc;
  REAL kappa;
  int print_aij, print_aijkl, print_aijklmn;
  REAL Dr;
  REAL normgamma;
  int length;
} param_list_t;

void ode_adams_bash_2_solve(REAL *t, REAL *x, REAL h, int do_many,
                            param_list_t *param,param_list_t *param_d);
void ode_adams_bash_4_solve(REAL *t, REAL *x, REAL h,
                            void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param,param_list_t *param_d),
                            param_list_t *param);
void ode_rk_4_solve(REAL *t, REAL *x, REAL h,
                    void derivs(REAL t, REAL *x, REAL *diffx, param_list_t *param,param_list_t *param_d),
                    param_list_t *param);
void compute_psidot(REAL* psidot, REAL* psi, param_list_t *param,param_list_t *param_d, int do_adams_bash_2, int nr_times);
void compute_psidot_koch(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_dd(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_dd_2(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_vd(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_vl(REAL* psidot, REAL* psi, param_list_t *param);
void compute_psidot_ard(REAL* psidot, REAL* psi, param_list_t *param);
void tensor2(REAL *psi, REAL a[3][3], param_list_t *param);
void tensor4(REAL *psi, REAL a[3][3][3][3], param_list_t *param);
void tensor6(REAL *psi, REAL a[3][3][3][3][3][3], param_list_t *param);
void compute_psidot_rsc(REAL* psidot, REAL* psi, param_list_t *param);
void reverse_tensor2(REAL a[3][3], REAL *psi);
void reverse_tensor4(REAL a[3][3][3][3], REAL *psi);
void reverse_tensor6(REAL a[3][3][3][3][3][3], REAL *psi);
void diagonalize_sym(int n, REAL *A, REAL *eval, REAL *evec);
int param_bool(const char *p);
int param_int(const char *p);
REAL param_REAL(const char *p);
int param_choice(const char *p, ...);
void set_param_filename(const char *f);
void set_param_verbose_level(int v);
void done_with_param();
void param_ignore(const char *p);
