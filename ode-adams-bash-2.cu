#define NO_STDIO
#include "spherical.h"

static REAL *diffx;
static REAL *olddiffx;
static REAL *tempx;
static int first = 1;

/* x_vec = a_vec + b*c_vec */
__global__ void madd(REAL *x_vec, REAL *a_vec, REAL b, REAL *c_vec, int length) {
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  if (i<length)
    x_vec[i] = a_vec[i] + b*c_vec[i];
  __syncthreads();
}

/* x_vec += b*c_vec */
__global__ void maddto(REAL *x_vec, REAL b, REAL *c_vec, int length) {
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  if (i<length)
    x_vec[i] += b*c_vec[i];
  __syncthreads();
}

/* x_vec += b*c_vec+d*e_vec */
/*
__global__ void maddto_twice(REAL *x_vec, REAL b, REAL *c_vec, REAL d, REAL *e_vec, int length) {
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  if (i<length)
    x_vec[i] += b*c_vec[i] + d*e_vec[i];
  __syncthreads();
}
*/

void ode_adams_bash_2_solve(REAL *t, REAL *x, REAL h, int do_many,
                            param_list_t *param, param_list_t *param_d) {
  if (first) {
    first = 0;
    cudaMalloc((void**)&diffx,sizeof(REAL)*param->length);
    cudaMalloc((void**)&tempx,sizeof(REAL)*param->length);
    cudaMalloc((void**)&olddiffx,sizeof(REAL)*param->length);
/* Midpoint method */
//    derivs(*t,x,diffx,param,param_d);
    compute_psidot(diffx,x,param,param_d,0,1);
    cudaMemcpy(olddiffx,diffx,param->length*sizeof(REAL),cudaMemcpyDeviceToDevice);
    madd<<<param->length/64+1,64>>>(tempx,x,h/2,diffx, param->length);
//    derivs(*t+h/2,tempx,diffx,param,param_d);
    compute_psidot(diffx,tempx,param,param_d,0,1);
    maddto<<<param->length/64+1,64>>>(x,h,diffx, param->length);
    cudaMemcpy(diffx,olddiffx,param->length*sizeof(REAL),cudaMemcpyDeviceToDevice);
  } else {
/* Adams-Bashforth method of order 2 */
//    derivs(*t,x,diffx,param,param_d);
    if (do_many)
      compute_psidot(diffx,x,param,param_d,1/*do_adams_bash_2*/,param->print_every);
    else
      compute_psidot(diffx,x,param,param_d,1/*do_adams_bash_2*/,1);
//    maddto_twice<<<param->length/64+1,64>>>(x,3*h/2,diffx,-h/2,olddiffx, param->length);
  }
  if (do_many)
    *t += param->print_every*h;
  else
    *t += h;
}
