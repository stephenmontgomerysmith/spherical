/*
 * Created from psidot-ard.conf by expand-iterate.pl.
 */

#include "spherical.h"
static int first_mc_0 = 1;
static REAL *mult_constant_0;
#define mc_count_0 39
static void initialize_mc_0(param_list_t *param);

/* Converted into a threaded program by expand-thread.pl. */

#include "threads.h"

struct job_info_s {int lstart, lend; sem_t sem_start, sem_end;};
static struct {
  REAL* psidot;
  REAL* psi;
  REAL Dr[3][3];
} thread_pass_args_0;
static thread_return_t thread_function_0(void *arg);
static int first_thread_0 = 1;
static struct job_info_s *job_info_0;

void compute_psidot_ard(REAL* psidot, REAL* psi, param_list_t *param) {
  REAL a[3][3], a2[3][3], gamm2[3][3];
  int i1,i2,i3;
  REAL Dr[3][3];

  tensor2(psi,a);
  memset(a2,0,sizeof(a2));
  memset(gamm2,0,sizeof(gamm2));
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++) for (i3=0;i3<3;i3++) {
    a2[i1][i3] += a[i1][i2]*a[i2][i3];
    gamm2[i1][i3] += param->gamm[i1*3+i2]*param->gamm[i2*3+i3];
  }
  memset(Dr,0,sizeof(Dr));
  for (i1=0;i1<3;i1++)
    Dr[i1][i1] = param->b1*param->normgamma;
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++)
    Dr[i1][i2] += param->b2*param->normgamma*a[i1][i2] + param->b3*param->normgamma*a2[i1][i2]
                  + param->b4/2*param->gamm[i1*3+i2]
                  + param->b5/4/param->normgamma*gamm2[i1][i2];

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
    if (abs(m-2)<=l-2)
      mc[0] = ((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[0] = 0;
    if (abs(m)<=l-2)
      mc[1] = -((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[1] = 0;
    if (abs(m+2)<=l-2)
      mc[2] = ((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[2] = 0;
    if (abs(m-2)<=l)
      mc[3] = ((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[3] = 0;
    if (abs(m)<=l)
      mc[4] = (-4*pow(ll,3)-2*pow(ll,4)-3*pow(mm,2)+pow(ll,2)*(1-2*pow(mm,2))+ll*(3-2*pow(mm,2)))/(-6+8*ll+8*pow(ll,2));
    else
      mc[4] = 0;
    if (abs(m+2)<=l)
      mc[5] = ((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[5] = 0;
    if (abs(m-2)<=l+2)
      mc[6] = (ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[6] = 0;
    if (abs(m)<=l+2)
      mc[7] = -(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[7] = 0;
    if (abs(m+2)<=l+2)
      mc[8] = (ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[8] = 0;
    if (abs(m-2)<=l-2)
      mc[9] = -((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[9] = 0;
    if (abs(m+2)<=l-2)
      mc[10] = ((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[10] = 0;
    if (abs(m-2)<=l)
      mc[11] = ((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(6-8*ll-8*pow(ll,2));
    else
      mc[11] = 0;
    if (abs(m+2)<=l)
      mc[12] = ((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[12] = 0;
    if (abs(m-2)<=l+2)
      mc[13] = -(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[13] = 0;
    if (abs(m+2)<=l+2)
      mc[14] = (ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[14] = 0;
    if (abs(m-1)<=l-2)
      mc[15] = -(((-2+ll)*(1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)));
    else
      mc[15] = 0;
    if (abs(m+1)<=l-2)
      mc[16] = ((-2+ll)*(1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[16] = 0;
    if (abs(m-1)<=l)
      mc[17] = ((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[17] = 0;
    if (abs(m+1)<=l)
      mc[18] = ((3+2*ll+2*pow(ll,2))*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[18] = 0;
    if (abs(m-1)<=l+2)
      mc[19] = (ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[19] = 0;
    if (abs(m+1)<=l+2)
      mc[20] = -((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)));
    else
      mc[20] = 0;
    if (abs(m-2)<=l-2)
      mc[21] = -((-2+ll)*(1+ll)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[21] = 0;
    if (abs(m)<=l-2)
      mc[22] = -((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(2.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[22] = 0;
    if (abs(m+2)<=l-2)
      mc[23] = -((-2+ll)*(1+ll)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(4.*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[23] = 0;
    if (abs(m-2)<=l)
      mc[24] = -((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[24] = 0;
    if (abs(m)<=l)
      mc[25] = (-4*pow(ll,3)-2*pow(ll,4)-3*pow(mm,2)+pow(ll,2)*(1-2*pow(mm,2))+ll*(3-2*pow(mm,2)))/(-6+8*ll+8*pow(ll,2));
    else
      mc[25] = 0;
    if (abs(m+2)<=l)
      mc[26] = -((3+2*ll+2*pow(ll,2))*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*(-3+4*ll+4*pow(ll,2)));
    else
      mc[26] = 0;
    if (abs(m-2)<=l+2)
      mc[27] = -(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[27] = 0;
    if (abs(m)<=l+2)
      mc[28] = -(ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(2.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[28] = 0;
    if (abs(m+2)<=l+2)
      mc[29] = -(ll*(3+ll)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(4.*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[29] = 0;
    if (abs(m-1)<=l-2)
      mc[30] = ((-2+ll)*(1+ll)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[30] = 0;
    if (abs(m+1)<=l-2)
      mc[31] = ((-2+ll)*(1+ll)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[31] = 0;
    if (abs(m-1)<=l)
      mc[32] = ((3+2*ll+2*pow(ll,2))*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(6-8*ll-8*pow(ll,2));
    else
      mc[32] = 0;
    if (abs(m+1)<=l)
      mc[33] = ((3+2*ll+2*pow(ll,2))*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[33] = 0;
    if (abs(m-1)<=l+2)
      mc[34] = -((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)));
    else
      mc[34] = 0;
    if (abs(m+1)<=l+2)
      mc[35] = -((ll*(3+ll)*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)));
    else
      mc[35] = 0;
    if (abs(m)<=l-2)
      mc[36] = ((-2+ll)*(1+ll)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[36] = 0;
    if (abs(m)<=l)
      mc[37] = (-4*pow(ll,3)-2*pow(ll,4)+3*pow(mm,2)+2*ll*pow(mm,2)+2*pow(ll,2)*(-1+pow(mm,2)))/(-3+4*ll+4*pow(ll,2));
    else
      mc[37] = 0;
    if (abs(m)<=l+2)
      mc[38] = (ll*(3+ll)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[38] = 0;
  }
}

static thread_return_t thread_function_0(void *arg) {
  struct job_info_s* job_info = (struct job_info_s*)arg;
  while (1) {
    REAL* psidot;
    REAL* psi;
    REAL Dr[3][3];
    int lstart, lend;
    int l,m;
    /* Wait until thread is told to start. */
    sem_wait(job_info->sem_start);
    lstart = job_info->lstart;
    lend = job_info->lend;
    /* Copy data to thread. */
    memcpy(&(psidot),&(thread_pass_args_0.psidot),sizeof(psidot));
    memcpy(&(psi),&(thread_pass_args_0.psi),sizeof(psi));
    memcpy(&(Dr),&(thread_pass_args_0.Dr),sizeof(Dr));
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      REAL *mc = mult_constant_0 + mc_count_0*mc_ind(l,m);
      psidot[ind(l,m,0)] += Dr[0][0]*(mc[0]*psi[ind(l-2,m-2,0)]+mc[1]*psi[ind(l-2,m,0)]+mc[2]*psi[ind(l-2,m+2,0)]+mc[3]*psi[ind(l,m-2,0)]+mc[4]*psi[ind(l,m,0)]+mc[5]*psi[ind(l,m+2,0)]+mc[6]*psi[ind(l+2,m-2,0)]+mc[7]*psi[ind(l+2,m,0)]+mc[8]*psi[ind(l+2,m+2,0)])
                        + Dr[0][1]*((-mc[9])*psi[ind(l-2,m-2,1)]+(-mc[10])*psi[ind(l-2,m+2,1)]+(-mc[11])*psi[ind(l,m-2,1)]+(-mc[12])*psi[ind(l,m+2,1)]+(-mc[13])*psi[ind(l+2,m-2,1)]+(-mc[14])*psi[ind(l+2,m+2,1)])
                        + Dr[0][2]*(mc[15]*psi[ind(l-2,m-1,0)]+mc[16]*psi[ind(l-2,m+1,0)]+mc[17]*psi[ind(l,m-1,0)]+mc[18]*psi[ind(l,m+1,0)]+mc[19]*psi[ind(l+2,m-1,0)]+mc[20]*psi[ind(l+2,m+1,0)])
                        + Dr[1][1]*(mc[21]*psi[ind(l-2,m-2,0)]+mc[22]*psi[ind(l-2,m,0)]+mc[23]*psi[ind(l-2,m+2,0)]+mc[24]*psi[ind(l,m-2,0)]+mc[25]*psi[ind(l,m,0)]+mc[26]*psi[ind(l,m+2,0)]+mc[27]*psi[ind(l+2,m-2,0)]+mc[28]*psi[ind(l+2,m,0)]+mc[29]*psi[ind(l+2,m+2,0)])
                        + Dr[1][2]*((-mc[30])*psi[ind(l-2,m-1,1)]+(-mc[31])*psi[ind(l-2,m+1,1)]+(-mc[32])*psi[ind(l,m-1,1)]+(-mc[33])*psi[ind(l,m+1,1)]+(-mc[34])*psi[ind(l+2,m-1,1)]+(-mc[35])*psi[ind(l+2,m+1,1)])
                        + Dr[2][2]*(mc[36]*psi[ind(l-2,m,0)]+mc[37]*psi[ind(l,m,0)]+mc[38]*psi[ind(l+2,m,0)]);
      psidot[ind(l,m,1)] += Dr[0][0]*(mc[0]*psi[ind(l-2,m-2,1)]+mc[1]*psi[ind(l-2,m,1)]+mc[2]*psi[ind(l-2,m+2,1)]+mc[3]*psi[ind(l,m-2,1)]+mc[4]*psi[ind(l,m,1)]+mc[5]*psi[ind(l,m+2,1)]+mc[6]*psi[ind(l+2,m-2,1)]+mc[7]*psi[ind(l+2,m,1)]+mc[8]*psi[ind(l+2,m+2,1)])
                        + Dr[0][1]*(mc[9]*psi[ind(l-2,m-2,0)]+mc[10]*psi[ind(l-2,m+2,0)]+mc[11]*psi[ind(l,m-2,0)]+mc[12]*psi[ind(l,m+2,0)]+mc[13]*psi[ind(l+2,m-2,0)]+mc[14]*psi[ind(l+2,m+2,0)])
                        + Dr[0][2]*(mc[15]*psi[ind(l-2,m-1,1)]+mc[16]*psi[ind(l-2,m+1,1)]+mc[17]*psi[ind(l,m-1,1)]+mc[18]*psi[ind(l,m+1,1)]+mc[19]*psi[ind(l+2,m-1,1)]+mc[20]*psi[ind(l+2,m+1,1)])
                        + Dr[1][1]*(mc[21]*psi[ind(l-2,m-2,1)]+mc[22]*psi[ind(l-2,m,1)]+mc[23]*psi[ind(l-2,m+2,1)]+mc[24]*psi[ind(l,m-2,1)]+mc[25]*psi[ind(l,m,1)]+mc[26]*psi[ind(l,m+2,1)]+mc[27]*psi[ind(l+2,m-2,1)]+mc[28]*psi[ind(l+2,m,1)]+mc[29]*psi[ind(l+2,m+2,1)])
                        + Dr[1][2]*(mc[30]*psi[ind(l-2,m-1,0)]+mc[31]*psi[ind(l-2,m+1,0)]+mc[32]*psi[ind(l,m-1,0)]+mc[33]*psi[ind(l,m+1,0)]+mc[34]*psi[ind(l+2,m-1,0)]+mc[35]*psi[ind(l+2,m+1,0)])
                        + Dr[2][2]*(mc[36]*psi[ind(l-2,m,1)]+mc[37]*psi[ind(l,m,1)]+mc[38]*psi[ind(l+2,m,1)]);
    }
    /* Broadcast that thread is finished. */
    sem_post(job_info->sem_end);
  }
  return(0);
}
