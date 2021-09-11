#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double gamm[3][3] = {{0,1,0},{1,0,0},{0,0,0}};
double w[3] = {0,0,-1};
double lambda = 0.9;
double kappa = 0.5;

void derivs(double t, double x[3], double diffx[3]) {
  int i,j;
  double gamm_pp;

  gamm_pp = 0;
  for (i=0;i<3;i++) for (j=0;j<3;j++)
    gamm_pp += gamm[i][j]*x[i]*x[j];

  for (i=0;i<3;i++) {
    diffx[i] = 0;
    for (j=0;j<3;j++)
      diffx[i] += gamm[i][j]*x[j];
    diffx[i] -= gamm_pp*x[i];
  }
  for (i=0;i<3;i++)
    diffx[i] *= lambda*(1+kappa*pow(gamm_pp,2));

  diffx[0] += w[1]*x[2]-w[2]*x[1];
  diffx[1] += w[2]*x[0]-w[0]*x[2];
  diffx[2] += w[0]*x[1]-w[1]*x[0];

  for (i=0;i<3;i++)
    diffx[i] *= 0.5;
}

double rand_norm() {
  double u1, u2;

  u1 = arc4random()/(double)((1ll<<32)-1);
  u2 = arc4random()/(double)((1ll<<32)-1);
  return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

double k1[3],k2[3],k3[3],k4[3],tempx[3];

void ode_rk4_solve(double *t, double *x, double h,
               void derivs(double t, double *x, double *diffx))
{
  int i;

  derivs(*t,x,k1);
  for (i=0;i<3;i++) tempx[i] = x[i] + 0.5*h*k1[i];
  derivs(*t+0.5*h,tempx,k2);
  for (i=0;i<3;i++) tempx[i] = x[i] + 0.5*h*k2[i];
  derivs(*t+0.5*h,tempx,k3);
  for (i=0;i<3;i++) tempx[i] = x[i] + h*k3[i];
  derivs(*t+h,tempx,k4);
  for (i=0;i<3;i++)
    x[i] += h/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
  *t += h;
}

int main() {
  double r1,r2,r3,norm;
  double x[3];
  double t;
  int i,j;
  double a11[4001];
  double a22[4001];
  double a33[4001];

  memset(a11,0,sizeof(a11));
  memset(a22,0,sizeof(a22));
  memset(a33,0,sizeof(a33));
  for (i=0;i<100000;i++) {
    r1 = rand_norm();
    r2 = rand_norm();
    r3 = rand_norm();
    x[0] = r1/sqrt(r1*r1+r2*r2+r3*r3);
    x[1] = r2/sqrt(r1*r1+r2*r2+r3*r3);
    x[2] = r3/sqrt(r1*r1+r2*r2+r3*r3);
    t = 0;
    a11[0] += x[0]*x[0];
    a22[0] += x[1]*x[1];
    a33[0] += x[2]*x[2];
    for (j=1;j<=4000;j++) {
      ode_rk4_solve(&t,x,0.01,derivs);
      norm = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      x[0] /= norm;
      x[1] /= norm;
      x[2] /= norm;
      a11[j] += x[0]*x[0];
      a22[j] += x[1]*x[1];
      a33[j] += x[2]*x[2];
    }
  }
  for (j=0;j<=4000;j+=10)
    printf("%g %g %g %g\n",0.01*j,a11[j]/100000,a22[j]/100000,a33[j]/100000);
  exit(0);
}
