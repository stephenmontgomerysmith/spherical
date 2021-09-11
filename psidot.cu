/*
 * Created from psidot.conf by expand-iterate.pl.
 */

#define NO_STDIO
#include "spherical.h"

#define block_width 8

static int first_mc_0 = 1;
static REAL *mult_constant_0;
#define mc_count_0 45
static void initialize_mc_0(param_list_t *param);

/* Declarations of threading go here. */
static int first_thread_0 = 1;
REAL *mult_constant_0_d;
texture<float,1,cudaReadModeElementType> mult_constant_0_r;

#define mc_0(i) mult_constant_0[((i)*param->data_width/2+l/2)*param->data_width+m]

//#define mc_0_r(i) tex1Dfetch(mult_constant_0_r,((i)*data_width/2+l/2)*data_width+m)

#define mc_0_r(i) mult_constant_0[((i)*data_width/2+l/2)*data_width+m]

#define promote_to_float2(x) (*((float2*)(&(x))))

__global__ void cuda_thread_0(REAL* psidot_io, REAL* psi_in, param_list_t *param, REAL* mult_constant_0, int do_adams_bash_2) {
  __shared__ REAL lambda;
  __shared__ REAL Dr;
  __shared__ REAL gamm[9], w[3];
#define index(psi,ll,mm,c) psi[(((ll)-l_base+PADDING)/2*(block_width+2*PADDING)+(mm)-m_base+PADDING)*2+c]
  __shared__ REAL psi[(block_width+PADDING)*(block_width+2*PADDING)*2];
  int l = 2*(blockIdx.y*blockDim.y + threadIdx.y);
  int m = blockIdx.x*blockDim.x + threadIdx.x;
  int l_base = 2*(blockIdx.y*blockDim.y);
  int m_base = blockIdx.x*blockDim.x;
/*if (threadIdx.y&(block_width/2)) l_half = PADDING; else l_half = -PADDING;
  if (threadIdx.x&(block_width/2)) m_half = PADDING; else m_half = -PADDING; */
  int l_half = ((threadIdx.y&(block_width/2))-(block_width/4))/(block_width/4)*PADDING;
  int m_half = ((threadIdx.x&(block_width/2))-(block_width/4))/(block_width/4)*PADDING;
  REAL psidot[2];
  int job_nr = threadIdx.y*block_width+threadIdx.x;
  int data_width = param->data_width;
  REAL h = param->h;
#undef ind
#define ind(l,m,c) ind_macro(l,m,c,data_width)
  float2 olddiffx;

  if (m_base<=l_base+block_width) {
    if (job_nr==0)
      lambda = param->lambda;
    if (job_nr==1)
      Dr = param->Dr;
    if (0<=job_nr-2 && job_nr-2<9)
      gamm[job_nr-2] = param->gamm[job_nr-2];
    if (0<=job_nr-11 && job_nr-11<3)
      w[job_nr-11] = param->w[job_nr-11];
    promote_to_float2(index(psi,l,m,0)) = promote_to_float2(psi_in[ind(l,m,0)]);
    promote_to_float2(index(psi,l,m+m_half,0)) = promote_to_float2(psi_in[ind(l,m+m_half,0)]);
    promote_to_float2(index(psi,l+l_half,m,0)) = promote_to_float2(psi_in[ind(l+l_half,m,0)]);
    promote_to_float2(index(psi,l+l_half,m+m_half,0)) = promote_to_float2(psi_in[ind(l+l_half,m+m_half,0)]);
/*
 * Condon-Shortley phase:
 * Y_l^(-m) = (-1)^m conj(Y_l^m)
 */
    if (m<=PADDING) {
      index(psi,l,-m,0) = (m&1) ? -index(psi,l,m,0) : index(psi,l,m,0);
      index(psi,l,-m,1) = !(m&1) ? -index(psi,l,m,1) : index(psi,l,m,1);
    }
//    promote_to_float2(psidot[0]) = promote_to_float2(psidot_io[ind(l,m,0)]);
  }
  __syncthreads();
  if (l<=param->max_order && m<=l) {
      psidot[0] = w[0]*((-mc_0_r(0))*index(psi,l,m-1,1)+(-mc_0_r(1))*index(psi,l,m+1,1))
                        + w[1]*(mc_0_r(2)*index(psi,l,m-1,0)+mc_0_r(3)*index(psi,l,m+1,0))
                        + w[2]*((-mc_0_r(4))*index(psi,l,m,1))
                        + lambda*gamm[0*3+0]*(mc_0_r(5)*index(psi,l-2,m-2,0)+mc_0_r(6)*index(psi,l-2,m,0)+mc_0_r(7)*index(psi,l-2,m+2,0)+mc_0_r(8)*index(psi,l,m-2,0)+mc_0_r(9)*index(psi,l,m,0)+mc_0_r(10)*index(psi,l,m+2,0)+mc_0_r(11)*index(psi,l+2,m-2,0)+mc_0_r(12)*index(psi,l+2,m,0)+mc_0_r(13)*index(psi,l+2,m+2,0))
                        + lambda*gamm[0*3+1]*((-mc_0_r(14))*index(psi,l-2,m-2,1)+(-mc_0_r(15))*index(psi,l-2,m+2,1)+(-mc_0_r(16))*index(psi,l,m-2,1)+(-mc_0_r(17))*index(psi,l,m+2,1)+(-mc_0_r(18))*index(psi,l+2,m-2,1)+(-mc_0_r(19))*index(psi,l+2,m+2,1))
                        + lambda*gamm[0*3+2]*(mc_0_r(20)*index(psi,l-2,m-1,0)+mc_0_r(21)*index(psi,l-2,m+1,0)+mc_0_r(22)*index(psi,l,m-1,0)+mc_0_r(23)*index(psi,l,m+1,0)+mc_0_r(24)*index(psi,l+2,m-1,0)+mc_0_r(25)*index(psi,l+2,m+1,0))
                        + lambda*gamm[1*3+1]*(mc_0_r(26)*index(psi,l-2,m-2,0)+mc_0_r(27)*index(psi,l-2,m,0)+mc_0_r(28)*index(psi,l-2,m+2,0)+mc_0_r(29)*index(psi,l,m-2,0)+mc_0_r(30)*index(psi,l,m,0)+mc_0_r(31)*index(psi,l,m+2,0)+mc_0_r(32)*index(psi,l+2,m-2,0)+mc_0_r(33)*index(psi,l+2,m,0)+mc_0_r(34)*index(psi,l+2,m+2,0))
                        + lambda*gamm[1*3+2]*((-mc_0_r(35))*index(psi,l-2,m-1,1)+(-mc_0_r(36))*index(psi,l-2,m+1,1)+(-mc_0_r(37))*index(psi,l,m-1,1)+(-mc_0_r(38))*index(psi,l,m+1,1)+(-mc_0_r(39))*index(psi,l+2,m-1,1)+(-mc_0_r(40))*index(psi,l+2,m+1,1))
                        + lambda*gamm[2*3+2]*(mc_0_r(41)*index(psi,l-2,m,0)+mc_0_r(42)*index(psi,l,m,0)+mc_0_r(43)*index(psi,l+2,m,0))
                        + Dr*(mc_0_r(44)*index(psi,l,m,0));
      psidot[1] = w[0]*(mc_0_r(0)*index(psi,l,m-1,0)+mc_0_r(1)*index(psi,l,m+1,0))
                        + w[1]*(mc_0_r(2)*index(psi,l,m-1,1)+mc_0_r(3)*index(psi,l,m+1,1))
                        + w[2]*(mc_0_r(4)*index(psi,l,m,0))
                        + lambda*gamm[0*3+0]*(mc_0_r(5)*index(psi,l-2,m-2,1)+mc_0_r(6)*index(psi,l-2,m,1)+mc_0_r(7)*index(psi,l-2,m+2,1)+mc_0_r(8)*index(psi,l,m-2,1)+mc_0_r(9)*index(psi,l,m,1)+mc_0_r(10)*index(psi,l,m+2,1)+mc_0_r(11)*index(psi,l+2,m-2,1)+mc_0_r(12)*index(psi,l+2,m,1)+mc_0_r(13)*index(psi,l+2,m+2,1))
                        + lambda*gamm[0*3+1]*(mc_0_r(14)*index(psi,l-2,m-2,0)+mc_0_r(15)*index(psi,l-2,m+2,0)+mc_0_r(16)*index(psi,l,m-2,0)+mc_0_r(17)*index(psi,l,m+2,0)+mc_0_r(18)*index(psi,l+2,m-2,0)+mc_0_r(19)*index(psi,l+2,m+2,0))
                        + lambda*gamm[0*3+2]*(mc_0_r(20)*index(psi,l-2,m-1,1)+mc_0_r(21)*index(psi,l-2,m+1,1)+mc_0_r(22)*index(psi,l,m-1,1)+mc_0_r(23)*index(psi,l,m+1,1)+mc_0_r(24)*index(psi,l+2,m-1,1)+mc_0_r(25)*index(psi,l+2,m+1,1))
                        + lambda*gamm[1*3+1]*(mc_0_r(26)*index(psi,l-2,m-2,1)+mc_0_r(27)*index(psi,l-2,m,1)+mc_0_r(28)*index(psi,l-2,m+2,1)+mc_0_r(29)*index(psi,l,m-2,1)+mc_0_r(30)*index(psi,l,m,1)+mc_0_r(31)*index(psi,l,m+2,1)+mc_0_r(32)*index(psi,l+2,m-2,1)+mc_0_r(33)*index(psi,l+2,m,1)+mc_0_r(34)*index(psi,l+2,m+2,1))
                        + lambda*gamm[1*3+2]*(mc_0_r(35)*index(psi,l-2,m-1,0)+mc_0_r(36)*index(psi,l-2,m+1,0)+mc_0_r(37)*index(psi,l,m-1,0)+mc_0_r(38)*index(psi,l,m+1,0)+mc_0_r(39)*index(psi,l+2,m-1,0)+mc_0_r(40)*index(psi,l+2,m+1,0))
                        + lambda*gamm[2*3+2]*(mc_0_r(41)*index(psi,l-2,m,1)+mc_0_r(42)*index(psi,l,m,1)+mc_0_r(43)*index(psi,l+2,m,1))
                        + Dr*(mc_0_r(44)*index(psi,l,m,1));
  }
  __syncthreads();
  if (l<=param->max_order && m<=l) {
    if (do_adams_bash_2) {
      olddiffx = promote_to_float2(psidot_io[ind(l,m,0)]);
      index(psi,l,m,0) += 3*h/2*psidot[0] - h/2*olddiffx.x;
      index(psi,l,m,1) += 3*h/2*psidot[1] - h/2*olddiffx.y;
      olddiffx =  promote_to_float2(psidot[0]);
      promote_to_float2(psi_in[ind(l,m,0)]) = promote_to_float2(index(psi,l,m,0));
    }
    promote_to_float2(psidot_io[ind(l,m,0)]) = promote_to_float2(psidot[0]);
  }
  __syncthreads();
}

void compute_psidot(REAL* psidot, REAL* psi, param_list_t *param, param_list_t *param_d, int do_adams_bash_2, int nr_times) {
  dim3 dimblock(block_width,block_width);
  dim3 dimgrid(param->data_width/block_width,param->data_width/block_width/2);
  int i;
  {
    if (first_mc_0) initialize_mc_0(param);
    if (first_thread_0) {
      first_thread_0 = 0;
      cudaMalloc((void**)&mult_constant_0_d, sizeof(REAL)*param->data_width/2*param->data_width*mc_count_0);
      cudaMemcpy(mult_constant_0_d,mult_constant_0,sizeof(REAL)*param->data_width/2*param->data_width*mc_count_0,cudaMemcpyHostToDevice);
      cudaBindTexture(0,mult_constant_0_r,mult_constant_0_d,sizeof(REAL)*param->data_width/2*param->data_width*mc_count_0);
    }

    for (i=0;i<nr_times;i++)
      cuda_thread_0<<<dimgrid,dimblock>>>(psidot, psi, param_d, mult_constant_0_d, do_adams_bash_2);
  }
}

static void initialize_mc_0(param_list_t *param) {
  int l,m;
  REAL ll,mm;

  first_mc_0 = 0;
  mult_constant_0 = (REAL*)malloc(sizeof(REAL)*param->data_width/2*param->data_width*mc_count_0);
  for (l=0;l<=param->max_order;l+=2) for (m=0;m<=l;m++) {
    ll = l;
    mm = m;
    if (abs(m-1)<=l)
      mc_0(0) = -(sqrt(1+ll-mm)*sqrt(ll+mm))/4.;
    else
      mc_0(0) = 0;
    if (abs(m+1)<=l)
      mc_0(1) = -(sqrt(ll-mm)*sqrt(1+ll+mm))/4.;
    else
      mc_0(1) = 0;
    if (abs(m-1)<=l)
      mc_0(2) = -(sqrt(1+ll-mm)*sqrt(ll+mm))/4.;
    else
      mc_0(2) = 0;
    if (abs(m+1)<=l)
      mc_0(3) = (sqrt(ll-mm)*sqrt(1+ll+mm))/4.;
    else
      mc_0(3) = 0;
    if (abs(m)<=l)
      mc_0(4) = -mm/2.;
    else
      mc_0(4) = 0;
    if (abs(m-2)<=l-2)
      mc_0(5) = ((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(5) = 0;
    if (abs(m)<=l-2)
      mc_0(6) = -((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(6) = 0;
    if (abs(m+2)<=l-2)
      mc_0(7) = ((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(7) = 0;
    if (abs(m-2)<=l)
      mc_0(8) = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(8) = 0;
    if (abs(m)<=l)
      mc_0(9) = (ll+pow(ll,2)-3*pow(mm,2))/(12-16*ll-16*pow(ll,2));
    else
      mc_0(9) = 0;
    if (abs(m+2)<=l)
      mc_0(10) = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(8.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(10) = 0;
    if (abs(m-2)<=l+2)
      mc_0(11) = -((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll)));
    else
      mc_0(11) = 0;
    if (abs(m)<=l+2)
      mc_0(12) = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc_0(12) = 0;
    if (abs(m+2)<=l+2)
      mc_0(13) = -((ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll)));
    else
      mc_0(13) = 0;
    if (abs(m-2)<=l-2)
      mc_0(14) = -((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(14) = 0;
    if (abs(m+2)<=l-2)
      mc_0(15) = ((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(15) = 0;
    if (abs(m-2)<=l)
      mc_0(16) = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(16) = 0;
    if (abs(m+2)<=l)
      mc_0(17) = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(17) = 0;
    if (abs(m-2)<=l+2)
      mc_0(18) = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc_0(18) = 0;
    if (abs(m+2)<=l+2)
      mc_0(19) = -((ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)));
    else
      mc_0(19) = 0;
    if (abs(m-1)<=l-2)
      mc_0(20) = -((1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(20) = 0;
    if (abs(m+1)<=l-2)
      mc_0(21) = ((1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(21) = 0;
    if (abs(m-1)<=l)
      mc_0(22) = (3*(1-2*mm)*sqrt(1+ll-mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(22) = 0;
    if (abs(m+1)<=l)
      mc_0(23) = (-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(23) = 0;
    if (abs(m-1)<=l+2)
      mc_0(24) = -((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc_0(24) = 0;
    if (abs(m+1)<=l+2)
      mc_0(25) = (ll*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll));
    else
      mc_0(25) = 0;
    if (abs(m-2)<=l-2)
      mc_0(26) = -((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(26) = 0;
    if (abs(m)<=l-2)
      mc_0(27) = -((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(27) = 0;
    if (abs(m+2)<=l-2)
      mc_0(28) = -((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(28) = 0;
    if (abs(m-2)<=l)
      mc_0(29) = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(29) = 0;
    if (abs(m)<=l)
      mc_0(30) = (ll+pow(ll,2)-3*pow(mm,2))/(12-16*ll-16*pow(ll,2));
    else
      mc_0(30) = 0;
    if (abs(m+2)<=l)
      mc_0(31) = (3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(8.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(31) = 0;
    if (abs(m-2)<=l+2)
      mc_0(32) = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll));
    else
      mc_0(32) = 0;
    if (abs(m)<=l+2)
      mc_0(33) = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc_0(33) = 0;
    if (abs(m+2)<=l+2)
      mc_0(34) = (ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll));
    else
      mc_0(34) = 0;
    if (abs(m-1)<=l-2)
      mc_0(35) = ((1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(35) = 0;
    if (abs(m+1)<=l-2)
      mc_0(36) = ((1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(36) = 0;
    if (abs(m-1)<=l)
      mc_0(37) = (3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(37) = 0;
    if (abs(m+1)<=l)
      mc_0(38) = (-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc_0(38) = 0;
    if (abs(m-1)<=l+2)
      mc_0(39) = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll));
    else
      mc_0(39) = 0;
    if (abs(m+1)<=l+2)
      mc_0(40) = (ll*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll));
    else
      mc_0(40) = 0;
    if (abs(m)<=l-2)
      mc_0(41) = ((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc_0(41) = 0;
    if (abs(m)<=l)
      mc_0(42) = (ll+pow(ll,2)-3*pow(mm,2))/(-6+8*ll+8*pow(ll,2));
    else
      mc_0(42) = 0;
    if (abs(m)<=l+2)
      mc_0(43) = -((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc_0(43) = 0;
    if (abs(m)<=l)
      mc_0(44) = -(ll*(1+ll));
    else
      mc_0(44) = 0;
  }
}
