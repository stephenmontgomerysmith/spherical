/*
 * Created from psidot-vl.conf by expand-iterate.pl.
 */

#include "spherical.h"
static int first_mc_0 = 1;
static REAL *mult_constant_0;
#define mc_count_0 45
static void initialize_mc_0(param_list_t *param);

/* Converted into a threaded program by expand-thread.pl. */

#include "threads.h"

struct job_info_s {int lstart, lend; sem_t sem_start, sem_end;};
static struct {
  REAL* psidot;
  REAL* psi;
  param_list_t *param;
  REAL s;
  REAL lambda;
  REAL Dr;
} thread_pass_args_0;
static thread_return_t thread_function_0(void *arg);
static int first_thread_0 = 1;
static struct job_info_s *job_info_0;

void compute_psidot_vl(REAL* psidot, REAL* psi, param_list_t *param) {
  REAL s;
  REAL a2[3][3];
  REAL lambda;
  int i1,i2;
  REAL Dr;

  Dr = param->CI * param->normgamma;

  s = 0;
  tensor2(psi,a2);
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++)
    s += a2[i1][i2]*param->gamm[i1*3+i2];
  lambda = evaluate_string(param->vl_fun,'s',s/param->normgamma);

  if (first_mc_0) initialize_mc_0(param);
  {
    int i;
    /* Start threads the first time around. */
    if (first_thread_0) {
      first_thread_0 = 0;
      job_info_0 = (struct job_info_s*)malloc(param->nr_threads*sizeof(struct job_info_s));
      for (i=0;i<param->nr_threads;i++) {
        sem_init(job_info_0[i].sem_start);
        sem_init(job_info_0[i].sem_end);
        /* The ranges for l are assigned to each thread.  Note that the work
           required to calculate for l between lstart and lend is
           proportional to lend^2 - lstart^2. */
        if (i==0)
          job_info_0[i].lstart = 0;
        else
          job_info_0[i].lstart = job_info_0[i-1].lend;
        if (i==param->nr_threads-1)
          job_info_0[i].lend=param->max_order+1;
        else
          job_info_0[i].lend = (REAL)(param->max_order+1)*sqrt((i+1.)/(REAL)param->nr_threads);
        job_info_0[i].lend += (job_info_0[i].lend)%2; /* round up to next even number. */
        thread_create(thread_function_0, job_info_0+i);
      }
    }
    /* Copy data to threads. */
    memcpy(&(thread_pass_args_0.psidot),&(psidot),sizeof(psidot));
    memcpy(&(thread_pass_args_0.psi),&(psi),sizeof(psi));
    memcpy(&(thread_pass_args_0.param),&(param),sizeof(param));
    memcpy(&(thread_pass_args_0.s),&(s),sizeof(s));
    memcpy(&(thread_pass_args_0.lambda),&(lambda),sizeof(lambda));
    memcpy(&(thread_pass_args_0.Dr),&(Dr),sizeof(Dr));
    /* Start the threads, and then wait for each thread to finish. */
    for (i=0;i<param->nr_threads;i++)
      sem_post(job_info_0[i].sem_start);
    for (i=0;i<param->nr_threads;i++)
      sem_wait(job_info_0[i].sem_end);
  }
}

static void initialize_mc_0(param_list_t *param) {
  int l,m;
  REAL ll,mm;
  REAL *mc;

  first_mc_0 = 0;
  mult_constant_0 = (REAL*)malloc(sizeof(REAL)*mc_ind(param->max_order+2,0)*mc_count_0);
  for (l=0;l<=param->max_order;l+=2) for (m=0;m<=l;m++) {
    ll = l;
    mm = m;
    mc = mult_constant_0 + mc_count_0*mc_ind(l,m);
    if (abs(m-1)<=l)
      mc[0] = -(sqrt(1+ll-mm)*sqrt(ll+mm))/4.;
    else
      mc[0] = 0;
    if (abs(m+1)<=l)
      mc[1] = -(sqrt(ll-mm)*sqrt(1+ll+mm))/4.;
    else
      mc[1] = 0;
    if (abs(m-1)<=l)
      mc[2] = -(sqrt(1+ll-mm)*sqrt(ll+mm))/4.;
    else
      mc[2] = 0;
    if (abs(m+1)<=l)
      mc[3] = (sqrt(ll-mm)*sqrt(1+ll+mm))/4.;
    else
      mc[3] = 0;
    if (abs(m)<=l)
      mc[4] = -mm/2.;
    else
      mc[4] = 0;
    if (abs(m-2)<=l-2)
      mc[5] = ((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[5] = 0;
    if (abs(m)<=l-2)
      mc[6] = -((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[6] = 0;
    if (abs(m+2)<=l-2)
      mc[7] = ((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[7] = 0;
    if (abs(m-2)<=l)
      mc[8] = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[8] = 0;
    if (abs(m)<=l)
      mc[9] = (ll+pow(ll,2)-3*pow(mm,2))/(12-16*ll-16*pow(ll,2));
    else
      mc[9] = 0;
    if (abs(m+2)<=l)
      mc[10] = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(8.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[10] = 0;
    if (abs(m-2)<=l+2)
      mc[11] = -((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll)));
    else
      mc[11] = 0;
    if (abs(m)<=l+2)
      mc[12] = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc[12] = 0;
    if (abs(m+2)<=l+2)
      mc[13] = -((ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll)));
    else
      mc[13] = 0;
    if (abs(m-2)<=l-2)
      mc[14] = -((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[14] = 0;
    if (abs(m+2)<=l-2)
      mc[15] = ((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[15] = 0;
    if (abs(m-2)<=l)
      mc[16] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[16] = 0;
    if (abs(m+2)<=l)
      mc[17] = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[17] = 0;
    if (abs(m-2)<=l+2)
      mc[18] = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc[18] = 0;
    if (abs(m+2)<=l+2)
      mc[19] = -((ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)));
    else
      mc[19] = 0;
    if (abs(m-1)<=l-2)
      mc[20] = -((1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[20] = 0;
    if (abs(m+1)<=l-2)
      mc[21] = ((1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[21] = 0;
    if (abs(m-1)<=l)
      mc[22] = (3*(1-2*mm)*sqrt(1+ll-mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[22] = 0;
    if (abs(m+1)<=l)
      mc[23] = (-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[23] = 0;
    if (abs(m-1)<=l+2)
      mc[24] = -((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc[24] = 0;
    if (abs(m+1)<=l+2)
      mc[25] = (ll*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll));
    else
      mc[25] = 0;
    if (abs(m-2)<=l-2)
      mc[26] = -((1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[26] = 0;
    if (abs(m)<=l-2)
      mc[27] = -((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[27] = 0;
    if (abs(m+2)<=l-2)
      mc[28] = -((1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(8.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[28] = 0;
    if (abs(m-2)<=l)
      mc[29] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[29] = 0;
    if (abs(m)<=l)
      mc[30] = (ll+pow(ll,2)-3*pow(mm,2))/(12-16*ll-16*pow(ll,2));
    else
      mc[30] = 0;
    if (abs(m+2)<=l)
      mc[31] = (3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(8.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[31] = 0;
    if (abs(m-2)<=l+2)
      mc[32] = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll));
    else
      mc[32] = 0;
    if (abs(m)<=l+2)
      mc[33] = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc[33] = 0;
    if (abs(m+2)<=l+2)
      mc[34] = (ll*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(24+16*ll));
    else
      mc[34] = 0;
    if (abs(m-1)<=l-2)
      mc[35] = ((1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[35] = 0;
    if (abs(m+1)<=l-2)
      mc[36] = ((1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[36] = 0;
    if (abs(m-1)<=l)
      mc[37] = (3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[37] = 0;
    if (abs(m+1)<=l)
      mc[38] = (-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[38] = 0;
    if (abs(m-1)<=l+2)
      mc[39] = (ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll));
    else
      mc[39] = 0;
    if (abs(m+1)<=l+2)
      mc[40] = (ll*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll));
    else
      mc[40] = 0;
    if (abs(m)<=l-2)
      mc[41] = ((1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[41] = 0;
    if (abs(m)<=l)
      mc[42] = (ll+pow(ll,2)-3*pow(mm,2))/(-6+8*ll+8*pow(ll,2));
    else
      mc[42] = 0;
    if (abs(m)<=l+2)
      mc[43] = -((ll*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc[43] = 0;
    if (abs(m)<=l)
      mc[44] = -(ll*(1+ll));
    else
      mc[44] = 0;
  }
}

static thread_return_t thread_function_0(void *arg) {
  struct job_info_s* job_info = (struct job_info_s*)arg;
  while (1) {
    REAL* psidot;
    REAL* psi;
    param_list_t *param;
    REAL s;
    REAL lambda;
    REAL Dr;
    int lstart, lend;
    int l,m;
    /* Wait until thread is told to start. */
    sem_wait(job_info->sem_start);
    lstart = job_info->lstart;
    lend = job_info->lend;
    /* Copy data to thread. */
    memcpy(&(psidot),&(thread_pass_args_0.psidot),sizeof(psidot));
    memcpy(&(psi),&(thread_pass_args_0.psi),sizeof(psi));
    memcpy(&(param),&(thread_pass_args_0.param),sizeof(param));
    memcpy(&(s),&(thread_pass_args_0.s),sizeof(s));
    memcpy(&(lambda),&(thread_pass_args_0.lambda),sizeof(lambda));
    memcpy(&(Dr),&(thread_pass_args_0.Dr),sizeof(Dr));
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      REAL *mc = mult_constant_0 + mc_count_0*mc_ind(l,m);
  /* Jeffery's Equation */
      psidot[ind(l,m,0)] += param->w[0]*((-mc[0])*psi[ind(l,m-1,1)]+(-mc[1])*psi[ind(l,m+1,1)])
                        + param->w[1]*(mc[2]*psi[ind(l,m-1,0)]+mc[3]*psi[ind(l,m+1,0)])
                        + param->w[2]*((-mc[4])*psi[ind(l,m,1)])
                        + lambda*param->gamm[0*3+0]*(mc[5]*psi[ind(l-2,m-2,0)]+mc[6]*psi[ind(l-2,m,0)]+mc[7]*psi[ind(l-2,m+2,0)]+mc[8]*psi[ind(l,m-2,0)]+mc[9]*psi[ind(l,m,0)]+mc[10]*psi[ind(l,m+2,0)]+mc[11]*psi[ind(l+2,m-2,0)]+mc[12]*psi[ind(l+2,m,0)]+mc[13]*psi[ind(l+2,m+2,0)])
                        + lambda*param->gamm[0*3+1]*((-mc[14])*psi[ind(l-2,m-2,1)]+(-mc[15])*psi[ind(l-2,m+2,1)]+(-mc[16])*psi[ind(l,m-2,1)]+(-mc[17])*psi[ind(l,m+2,1)]+(-mc[18])*psi[ind(l+2,m-2,1)]+(-mc[19])*psi[ind(l+2,m+2,1)])
                        + lambda*param->gamm[0*3+2]*(mc[20]*psi[ind(l-2,m-1,0)]+mc[21]*psi[ind(l-2,m+1,0)]+mc[22]*psi[ind(l,m-1,0)]+mc[23]*psi[ind(l,m+1,0)]+mc[24]*psi[ind(l+2,m-1,0)]+mc[25]*psi[ind(l+2,m+1,0)])
                        + lambda*param->gamm[1*3+1]*(mc[26]*psi[ind(l-2,m-2,0)]+mc[27]*psi[ind(l-2,m,0)]+mc[28]*psi[ind(l-2,m+2,0)]+mc[29]*psi[ind(l,m-2,0)]+mc[30]*psi[ind(l,m,0)]+mc[31]*psi[ind(l,m+2,0)]+mc[32]*psi[ind(l+2,m-2,0)]+mc[33]*psi[ind(l+2,m,0)]+mc[34]*psi[ind(l+2,m+2,0)])
                        + lambda*param->gamm[1*3+2]*((-mc[35])*psi[ind(l-2,m-1,1)]+(-mc[36])*psi[ind(l-2,m+1,1)]+(-mc[37])*psi[ind(l,m-1,1)]+(-mc[38])*psi[ind(l,m+1,1)]+(-mc[39])*psi[ind(l+2,m-1,1)]+(-mc[40])*psi[ind(l+2,m+1,1)])
                        + lambda*param->gamm[2*3+2]*(mc[41]*psi[ind(l-2,m,0)]+mc[42]*psi[ind(l,m,0)]+mc[43]*psi[ind(l+2,m,0)])
  /* Folgar-Tucker correction */
                        + Dr*(mc[44]*psi[ind(l,m,0)]);
  /* Jeffery's Equation */
      psidot[ind(l,m,1)] += param->w[0]*(mc[0]*psi[ind(l,m-1,0)]+mc[1]*psi[ind(l,m+1,0)])
                        + param->w[1]*(mc[2]*psi[ind(l,m-1,1)]+mc[3]*psi[ind(l,m+1,1)])
                        + param->w[2]*(mc[4]*psi[ind(l,m,0)])
                        + lambda*param->gamm[0*3+0]*(mc[5]*psi[ind(l-2,m-2,1)]+mc[6]*psi[ind(l-2,m,1)]+mc[7]*psi[ind(l-2,m+2,1)]+mc[8]*psi[ind(l,m-2,1)]+mc[9]*psi[ind(l,m,1)]+mc[10]*psi[ind(l,m+2,1)]+mc[11]*psi[ind(l+2,m-2,1)]+mc[12]*psi[ind(l+2,m,1)]+mc[13]*psi[ind(l+2,m+2,1)])
                        + lambda*param->gamm[0*3+1]*(mc[14]*psi[ind(l-2,m-2,0)]+mc[15]*psi[ind(l-2,m+2,0)]+mc[16]*psi[ind(l,m-2,0)]+mc[17]*psi[ind(l,m+2,0)]+mc[18]*psi[ind(l+2,m-2,0)]+mc[19]*psi[ind(l+2,m+2,0)])
                        + lambda*param->gamm[0*3+2]*(mc[20]*psi[ind(l-2,m-1,1)]+mc[21]*psi[ind(l-2,m+1,1)]+mc[22]*psi[ind(l,m-1,1)]+mc[23]*psi[ind(l,m+1,1)]+mc[24]*psi[ind(l+2,m-1,1)]+mc[25]*psi[ind(l+2,m+1,1)])
                        + lambda*param->gamm[1*3+1]*(mc[26]*psi[ind(l-2,m-2,1)]+mc[27]*psi[ind(l-2,m,1)]+mc[28]*psi[ind(l-2,m+2,1)]+mc[29]*psi[ind(l,m-2,1)]+mc[30]*psi[ind(l,m,1)]+mc[31]*psi[ind(l,m+2,1)]+mc[32]*psi[ind(l+2,m-2,1)]+mc[33]*psi[ind(l+2,m,1)]+mc[34]*psi[ind(l+2,m+2,1)])
                        + lambda*param->gamm[1*3+2]*(mc[35]*psi[ind(l-2,m-1,0)]+mc[36]*psi[ind(l-2,m+1,0)]+mc[37]*psi[ind(l,m-1,0)]+mc[38]*psi[ind(l,m+1,0)]+mc[39]*psi[ind(l+2,m-1,0)]+mc[40]*psi[ind(l+2,m+1,0)])
                        + lambda*param->gamm[2*3+2]*(mc[41]*psi[ind(l-2,m,1)]+mc[42]*psi[ind(l,m,1)]+mc[43]*psi[ind(l+2,m,1)])
  /* Folgar-Tucker correction */
                        + Dr*(mc[44]*psi[ind(l,m,1)]);
    }
    /* Broadcast that thread is finished. */
    sem_post(job_info->sem_end);
  }
  return(0);
}
