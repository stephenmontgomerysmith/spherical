/*
 * Created from psidot-vd.conf by expand-iterate.pl.
 */

#include "spherical.h"

extern int max_order;
static int first_mc_0 = 1;
static COMPLEX *mult_constant_0;
#define mc_count_0 1
static void initialize_mc_0();

/* Converted into a threaded program by expand-thread.pl. */
#include <pthread.h>
extern int nr_threads;
static pthread_mutex_t job_m = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t job_c = PTHREAD_COND_INITIALIZER;
struct job_info_s {int control, lstart, lend;};
static struct {
  COMPLEX* psidot;
  COMPLEX* psi;
  double Dr2;
} thread_pass_args_0;
static void *thread_function_0(void *arg);
static int first_thread_0 = 1;
static struct job_info_s *job_info_0;

extern double gamm[3][3];
extern double normgamma;
extern double C1;
extern int max_order;
extern int do_vd;

void compute_psidot_vd(COMPLEX* psidot, COMPLEX* psi) {
  int i1,i2,i3,i4;
  double a2[3][3];
  double a4[3][3][3][3];  double s;
  double Dr2;

  s = 0;
  if (do_vd == 1) {
    tensor2(psi,a2);
    for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++)
      s += a2[i1][i2]*gamm[i1][i2];
    s = fabs(s);
  } else if (do_vd ==2) {
    tensor4(psi,a4);
    for (i1=0;i1<3;i1++) for (i2=0;i2<3;i2++) for (i3=0;i3<3;i3++) for (i4=0;i4<3;i4++)
      s += a4[i1][i2][i3][i4]*gamm[i1][i2]*gamm[i3][i4];
    s = sqrt(s);
  }
  Dr2 = C1*s;

  {
    if (first_mc_0) initialize_mc_0();
    pthread_t pid;
    int i;
    /* Start threads the first time around. */
    if (first_thread_0) {
      first_thread_0 = 0;
      job_info_0 = malloc(nr_threads*sizeof(struct job_info_s));
      for (i=0;i<nr_threads;i++) {
        job_info_0[i].control = 0;
        /* The ranges for l are assigned to each thread.  Note that the work
           required to calculate for l between lstart and lend is
           proportional to lend^2 - lstart^2. */
        if (i==0)
          job_info_0[i].lstart = 0;
        else
          job_info_0[i].lstart = job_info_0[i-1].lend;
        if (i==nr_threads-1)
          job_info_0[i].lend=max_order+1;
        else
          job_info_0[i].lend = (double)(max_order+1)*sqrt((i+1.)/(double)nr_threads);
        job_info_0[i].lend += (job_info_0[i].lend)%2; /* round up to next even number. */
        pthread_create(&pid, NULL, thread_function_0, job_info_0+i);
        pthread_detach(pid);
      }
    }
    /* Copy data to threads. */
    memcpy(&(thread_pass_args_0.psidot),&(psidot),sizeof(psidot));
    memcpy(&(thread_pass_args_0.psi),&(psi),sizeof(psi));
    memcpy(&(thread_pass_args_0.Dr2),&(Dr2),sizeof(Dr2));
    /* Start each thread, and then wait for each thread to finish. */
    pthread_mutex_lock(&job_m);
    for (i=0;i<nr_threads;i++)
      job_info_0[i].control = 1;
    pthread_mutex_unlock(&job_m);
    pthread_cond_broadcast(&job_c);
    pthread_mutex_lock(&job_m);
    for (i=0;i<nr_threads;i++)
      while (job_info_0[i].control != 0)
        pthread_cond_wait(&job_c, &job_m);
    pthread_mutex_unlock(&job_m);
  }
}

static void initialize_mc_0() {
  int l,m;
  double ll,mm;
  COMPLEX *mc;

  first_mc_0 = 0;
  mult_constant_0 = malloc(sizeof(COMPLEX)*length*mc_count_0);
  for (l=0;l<=max_order;l+=2) for (m=0;m<=l;m++) {
    ll = l;
    mm = m;
    mc = mult_constant_0 + mc_count_0*ind(l,m);
    if (abs(m)<=l)
      mc[0] = (-(ll*(1+ll)))+(0)*I;
    else
      mc[0] = 0;
  }
}

static void * thread_function_0(void *arg) {
  struct job_info_s* job_info = arg;
  while (1) {
    int lstart, lend;
    /* Wait until thread is told to start. */
    pthread_mutex_lock(&job_m);
    while (job_info->control != 1)
      pthread_cond_wait(&job_c, &job_m);
    lstart = job_info->lstart;
    lend = job_info->lend;
    pthread_mutex_unlock(&job_m);
    COMPLEX* psidot;
    COMPLEX* psi;
    double Dr2;
    /* Copy data to thread. */
    memcpy(&(psidot),&(thread_pass_args_0.psidot),sizeof(psidot));
    memcpy(&(psi),&(thread_pass_args_0.psi),sizeof(psi));
    memcpy(&(Dr2),&(thread_pass_args_0.Dr2),sizeof(Dr2));
    int l,m;
    COMPLEX *mc;
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      mc = mult_constant_0 + mc_count_0*ind(l,m);
      psidot[ind(l,m)] += Dr2*(mc[0]*index(psi,l,m));
    }
    /* Broadcast that thread is finished. */
    pthread_mutex_lock(&job_m);
    job_info->control = 0;
    pthread_mutex_unlock(&job_m);
    pthread_cond_broadcast(&job_c);
  }
  return(NULL);
}
