/*
 * Created from psidot-vd.conf by expand-iterate.pl.
 */

#include "spherical.h"
static int first_mc_0 = 1;
static REAL *mult_constant_0;
#define mc_count_0 1
static void initialize_mc_0(param_list_t *param);

/* Converted into a threaded program by expand-thread.pl. */

#include "threads.h"

struct job_info_s {int lstart, lend; sem_t sem_start, sem_end;};
static struct {
  REAL* psidot;
  REAL* psi;
  REAL Dr2;
} thread_pass_args_0;
static thread_return_t thread_function_0(void *arg);
static int first_thread_0 = 1;
static struct job_info_s *job_info_0;

void compute_psidot_vd(REAL* psidot, REAL* psi, param_list_t *param) {
  int i1,i2;
  REAL a2[3][3];
  REAL s;
  REAL Dr2;

  s = 0;
  tensor2(psi,a2);
  for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++)
    s += a2[i1][i2]*param->gamm[i1*3+i2];
  Dr2 = param->normgamma*evaluate_string(param->vd_fun,'s',s/param->normgamma);

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
    memcpy(&(thread_pass_args_0.Dr2),&(Dr2),sizeof(Dr2));
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
    if (abs(m)<=l)
      mc[0] = -(ll*(1+ll));
    else
      mc[0] = 0;
  }
}

static thread_return_t thread_function_0(void *arg) {
  struct job_info_s* job_info = (struct job_info_s*)arg;
  while (1) {
    REAL* psidot;
    REAL* psi;
    REAL Dr2;
    int lstart, lend;
    int l,m;
    /* Wait until thread is told to start. */
    sem_wait(job_info->sem_start);
    lstart = job_info->lstart;
    lend = job_info->lend;
    /* Copy data to thread. */
    memcpy(&(psidot),&(thread_pass_args_0.psidot),sizeof(psidot));
    memcpy(&(psi),&(thread_pass_args_0.psi),sizeof(psi));
    memcpy(&(Dr2),&(thread_pass_args_0.Dr2),sizeof(Dr2));
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      REAL *mc = mult_constant_0 + mc_count_0*mc_ind(l,m);
      psidot[ind(l,m,0)] += Dr2*(mc[0]*psi[ind(l,m,0)]);
      psidot[ind(l,m,1)] += Dr2*(mc[0]*psi[ind(l,m,1)]);
    }
    /* Broadcast that thread is finished. */
    sem_post(job_info->sem_end);
  }
  return(0);
}
