#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

#ifdef USE_COMPLEX_DOUBLE
#define COMPLEX complex double
#else
#define COMPLEX complex
#endif

/*
 * psi[ind(l,m)] represents the Y_l^m spherical harmonic coefficient of psi.
 */
#define ind(l,m) ((l)*(l)/4+(m))

/* 
 * Y_l^(-m) = (-1)^m conj(Y_l^m)
 */
#define index(v,l,m) (                      \
  ((l)>max_order) ? 0 :                     \
  (abs(m)>(l))    ? 0 :                     \
  ((m)>=0)        ? v[ind(l,m)] :           \
  ((m)%2)         ? -conj(v[ind(l,-(m))]) : \
                    conj(v[ind(l,-(m))]) )

#define length ind(max_order+2,0)

void ode_adams_bash_2_solve(double *t, COMPLEX *x, double h,
                            void derivs(double t, COMPLEX *x, COMPLEX *diffx));
void ode_adams_bash_4_solve(double *t, COMPLEX *x, double h,
                            void derivs(double t, COMPLEX *x, COMPLEX *diffx));
void ode_rk_4_solve(double *t, COMPLEX *x, double h,
                    void derivs(double t, COMPLEX *x, COMPLEX *diffx));
void compute_psidot(COMPLEX* psidot, COMPLEX* psi);
void compute_psidot_koch(COMPLEX* psidot, COMPLEX* psi);
void compute_psidot_dd(COMPLEX* psidot, COMPLEX* psi);
void compute_psidot_dd_2(COMPLEX* psidot, COMPLEX* psi);
void compute_psidot_vd(COMPLEX* psidot, COMPLEX* psi);
void compute_psidot_vl(COMPLEX* psidot, COMPLEX* psi);
void compute_psidot_ard(COMPLEX* psidot, COMPLEX* psi);
void tensor2(COMPLEX *psi, double a[3][3]);
void tensor4(COMPLEX *psi, double a[3][3][3][3]);
void tensor6(COMPLEX *psi, double a[3][3][3][3][3][3]);
void compute_psidot_rsc(COMPLEX* psidot, COMPLEX* psi);
void reverse_tensor2(double a[3][3], COMPLEX *psi);
void reverse_tensor4(double a[3][3][3][3], COMPLEX *psi);
void reverse_tensor6(double a[3][3][3][3][3][3], COMPLEX *psi);
void diagonalize_sym(int n, double A[n][n], double eval[n], double evec[n][n]);
int param_bool(char *p);
int param_int(char *p);
double param_double(char *p);
int param_choice(char *p, ...);
void set_param_filename(char *f);
void set_param_verbose_level(int v);
void done_with_param();
void param_ignore(char *p);
