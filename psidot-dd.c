/*
 * Created from psidot-dd.conf by expand-iterate.pl.
 */

#include "spherical.h"
static int first_mc_0 = 1;
static REAL *mult_constant_0;
#define mc_count_0 289
static void initialize_mc_0(param_list_t *param);

/* Converted into a threaded program by expand-thread.pl. */

#include "threads.h"

struct job_info_s {int lstart, lend; sem_t sem_start, sem_end;};
static struct {
  REAL* psidot;
  REAL* psi;
  param_list_t *param;
  REAL a2[3][3];
  REAL a4[3][3][3][3];
} thread_pass_args_0;
static thread_return_t thread_function_0(void *arg);
static int first_thread_0 = 1;
static struct job_info_s *job_info_0;

void compute_psidot_dd(REAL* psidot, REAL* psi, param_list_t *param) {
  REAL a2[3][3], a4[3][3][3][3];

  tensor2(psi,a2);
  tensor4(psi,a4);

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
    memcpy(&(thread_pass_args_0.a2),&(a2),sizeof(a2));
    memcpy(&(thread_pass_args_0.a4),&(a4),sizeof(a4));
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
      mc[0] = -((sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)));
    else
      mc[0] = 0;
    if (abs(m)<=l-2)
      mc[1] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll));
    else
      mc[1] = 0;
    if (abs(m+2)<=l-2)
      mc[2] = -((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)));
    else
      mc[2] = 0;
    if (abs(m-2)<=l)
      mc[3] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(6-8*ll-8*pow(ll,2));
    else
      mc[3] = 0;
    if (abs(m)<=l)
      mc[4] = (-1+ll+pow(ll,2)+pow(mm,2))/(-3+4*ll+4*pow(ll,2));
    else
      mc[4] = 0;
    if (abs(m+2)<=l)
      mc[5] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(6-8*ll-8*pow(ll,2));
    else
      mc[5] = 0;
    if (abs(m-2)<=l+2)
      mc[6] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc[6] = 0;
    if (abs(m)<=l+2)
      mc[7] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc[7] = 0;
    if (abs(m+2)<=l+2)
      mc[8] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc[8] = 0;
    if (abs(m-2)<=l-2)
      mc[9] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll));
    else
      mc[9] = 0;
    if (abs(m+2)<=l-2)
      mc[10] = -((sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)));
    else
      mc[10] = 0;
    if (abs(m-2)<=l)
      mc[11] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[11] = 0;
    if (abs(m+2)<=l)
      mc[12] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(6-8*ll-8*pow(ll,2));
    else
      mc[12] = 0;
    if (abs(m-2)<=l+2)
      mc[13] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)));
    else
      mc[13] = 0;
    if (abs(m+2)<=l+2)
      mc[14] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll));
    else
      mc[14] = 0;
    if (abs(m-1)<=l-2)
      mc[15] = (sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll));
    else
      mc[15] = 0;
    if (abs(m+1)<=l-2)
      mc[16] = -((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)));
    else
      mc[16] = 0;
    if (abs(m-1)<=l)
      mc[17] = ((1-2*mm)*sqrt(1+ll-mm)*sqrt(ll+mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[17] = 0;
    if (abs(m+1)<=l)
      mc[18] = (sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(6-8*ll-8*pow(ll,2));
    else
      mc[18] = 0;
    if (abs(m-1)<=l+2)
      mc[19] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll));
    else
      mc[19] = 0;
    if (abs(m+1)<=l+2)
      mc[20] = -((sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc[20] = 0;
    if (abs(m-2)<=l-2)
      mc[21] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll));
    else
      mc[21] = 0;
    if (abs(m)<=l-2)
      mc[22] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll));
    else
      mc[22] = 0;
    if (abs(m+2)<=l-2)
      mc[23] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/((4-8*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll));
    else
      mc[23] = 0;
    if (abs(m-2)<=l)
      mc[24] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[24] = 0;
    if (abs(m)<=l)
      mc[25] = (-1+ll+pow(ll,2)+pow(mm,2))/(-3+4*ll+4*pow(ll,2));
    else
      mc[25] = 0;
    if (abs(m+2)<=l)
      mc[26] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[26] = 0;
    if (abs(m-2)<=l+2)
      mc[27] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)));
    else
      mc[27] = 0;
    if (abs(m)<=l+2)
      mc[28] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc[28] = 0;
    if (abs(m+2)<=l+2)
      mc[29] = -((sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(12+8*ll)));
    else
      mc[29] = 0;
    if (abs(m-1)<=l-2)
      mc[30] = -((sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)));
    else
      mc[30] = 0;
    if (abs(m+1)<=l-2)
      mc[31] = -((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/((2-4*ll)*sqrt(-3+2*ll)*sqrt(1+2*ll)));
    else
      mc[31] = 0;
    if (abs(m-1)<=l)
      mc[32] = (sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm))/(-6+8*ll+8*pow(ll,2));
    else
      mc[32] = 0;
    if (abs(m+1)<=l)
      mc[33] = (sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm))/(6-8*ll-8*pow(ll,2));
    else
      mc[33] = 0;
    if (abs(m-1)<=l+2)
      mc[34] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc[34] = 0;
    if (abs(m+1)<=l+2)
      mc[35] = -((sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(5+2*ll)*(6+4*ll)));
    else
      mc[35] = 0;
    if (abs(m)<=l-2)
      mc[36] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll));
    else
      mc[36] = 0;
    if (abs(m)<=l)
      mc[37] = (-1+2*ll+2*pow(ll,2)-2*pow(mm,2))/(-3+4*ll+4*pow(ll,2));
    else
      mc[37] = 0;
    if (abs(m)<=l+2)
      mc[38] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll));
    else
      mc[38] = 0;
    if (abs(m-2)<=l-2)
      mc[39] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-5+2*pow(ll,2)-3*mm-2*ll*mm+2*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[39] = 0;
    if (abs(m-4)<=l-2)
      mc[40] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3)));
    else
      mc[40] = 0;
    if (abs(m)<=l-2)
      mc[41] = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-4-ll+pow(ll,2)+pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[41] = 0;
    if (abs(m+2)<=l-2)
      mc[42] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(-5+2*pow(ll,2)+3*mm+2*ll*mm+2*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[42] = 0;
    if (abs(m+4)<=l-2)
      mc[43] = -(sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3)));
    else
      mc[43] = 0;
    if (abs(m-2)<=l-4)
      mc[44] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3)));
    else
      mc[44] = 0;
    if (abs(m-4)<=l-4)
      mc[45] = (sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[45] = 0;
    if (abs(m)<=l-4)
      mc[46] = (3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[46] = 0;
    if (abs(m+2)<=l-4)
      mc[47] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3)));
    else
      mc[47] = 0;
    if (abs(m+4)<=l-4)
      mc[48] = (sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[48] = 0;
    if (abs(m-2)<=l)
      mc[49] = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+ll+pow(ll,2)-2*mm+pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[49] = 0;
    if (abs(m-4)<=l)
      mc[50] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[50] = 0;
    if (abs(m)<=l)
      mc[51] = (18*pow(ll,3)+9*pow(ll,4)+6*ll*(-7+pow(mm,2))+pow(ll,2)*(-33+6*pow(mm,2))+9*(4-5*pow(mm,2)+pow(mm,4)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[51] = 0;
    if (abs(m+2)<=l)
      mc[52] = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+ll+pow(ll,2)+2*mm+pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[52] = 0;
    if (abs(m+4)<=l)
      mc[53] = (3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[53] = 0;
    if (abs(m-2)<=l+2)
      mc[54] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(-3+2*pow(ll,2)-mm+2*pow(mm,2)+2*ll*(2+mm)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[54] = 0;
    if (abs(m-4)<=l+2)
      mc[55] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[55] = 0;
    if (abs(m)<=l+2)
      mc[56] = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-2+3*ll+pow(ll,2)+pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[56] = 0;
    if (abs(m+2)<=l+2)
      mc[57] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-3+2*pow(ll,2)-2*ll*(-2+mm)+mm+2*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[57] = 0;
    if (abs(m+4)<=l+2)
      mc[58] = -(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[58] = 0;
    if (abs(m-2)<=l+4)
      mc[59] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))));
    else
      mc[59] = 0;
    if (abs(m-4)<=l+4)
      mc[60] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[60] = 0;
    if (abs(m)<=l+4)
      mc[61] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[61] = 0;
    if (abs(m+2)<=l+4)
      mc[62] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))));
    else
      mc[62] = 0;
    if (abs(m+4)<=l+4)
      mc[63] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[63] = 0;
    if (abs(m-2)<=l-2)
      mc[64] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-2*pow(ll,2)+3*mm+2*ll*mm-2*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[64] = 0;
    if (abs(m-4)<=l-2)
      mc[65] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*sqrt(1+2*ll)*(60-104*ll-48*pow(ll,2)+32*pow(ll,3)));
    else
      mc[65] = 0;
    if (abs(m+2)<=l-2)
      mc[66] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(-5+2*pow(ll,2)+3*mm+2*ll*mm+2*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[66] = 0;
    if (abs(m+4)<=l-2)
      mc[67] = -(sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3)));
    else
      mc[67] = 0;
    if (abs(m-2)<=l-4)
      mc[68] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[68] = 0;
    if (abs(m-4)<=l-4)
      mc[69] = (sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(240-736*ll+576*pow(ll,2)-128*pow(ll,3)));
    else
      mc[69] = 0;
    if (abs(m+2)<=l-4)
      mc[70] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[70] = 0;
    if (abs(m+4)<=l-4)
      mc[71] = (sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[71] = 0;
    if (abs(m-2)<=l)
      mc[72] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+ll+pow(ll,2)-2*mm+pow(mm,2)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[72] = 0;
    if (abs(m-4)<=l)
      mc[73] = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[73] = 0;
    if (abs(m+2)<=l)
      mc[74] = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+ll+pow(ll,2)+2*mm+pow(mm,2)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[74] = 0;
    if (abs(m+4)<=l)
      mc[75] = (3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[75] = 0;
    if (abs(m-2)<=l+2)
      mc[76] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(3-2*pow(ll,2)+mm-2*pow(mm,2)-2*ll*(2+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[76] = 0;
    if (abs(m-4)<=l+2)
      mc[77] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[77] = 0;
    if (abs(m+2)<=l+2)
      mc[78] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-3+2*pow(ll,2)-2*ll*(-2+mm)+mm+2*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[78] = 0;
    if (abs(m+4)<=l+2)
      mc[79] = -(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[79] = 0;
    if (abs(m-2)<=l+4)
      mc[80] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[80] = 0;
    if (abs(m-4)<=l+4)
      mc[81] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[81] = 0;
    if (abs(m+2)<=l+4)
      mc[82] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[82] = 0;
    if (abs(m+4)<=l+4)
      mc[83] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[83] = 0;
    if (abs(m-1)<=l-2)
      mc[84] = (3*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(7+ll-2*pow(ll,2)+3*mm+2*ll*mm-4*pow(mm,2)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[84] = 0;
    if (abs(m-3)<=l-2)
      mc[85] = ((5+2*ll-4*mm)*sqrt(1+ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[85] = 0;
    if (abs(m+1)<=l-2)
      mc[86] = (3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)+3*mm+4*pow(mm,2)+ll*(-1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[86] = 0;
    if (abs(m+3)<=l-2)
      mc[87] = -(sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*(5+2*ll+4*mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[87] = 0;
    if (abs(m-1)<=l-4)
      mc[88] = (3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[88] = 0;
    if (abs(m-3)<=l-4)
      mc[89] = (sqrt(ll-mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[89] = 0;
    if (abs(m+1)<=l-4)
      mc[90] = (-3*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[90] = 0;
    if (abs(m+3)<=l-4)
      mc[91] = (sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[91] = 0;
    if (abs(m-1)<=l)
      mc[92] = (-3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-6+ll+pow(ll,2)-3*mm+3*pow(mm,2)))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[92] = 0;
    if (abs(m-3)<=l)
      mc[93] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[93] = 0;
    if (abs(m+1)<=l)
      mc[94] = (-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(ll+pow(ll,2)+3*(-2+mm+pow(mm,2))))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[94] = 0;
    if (abs(m+3)<=l)
      mc[95] = (3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[95] = 0;
    if (abs(m-1)<=l+2)
      mc[96] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-4+2*pow(ll,2)-mm+4*pow(mm,2)+ll*(5+2*mm)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[96] = 0;
    if (abs(m-3)<=l+2)
      mc[97] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(ll+mm)*(-3+2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[97] = 0;
    if (abs(m+1)<=l+2)
      mc[98] = (-3*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(-4+2*pow(ll,2)+ll*(5-2*mm)+mm+4*pow(mm,2)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[98] = 0;
    if (abs(m+3)<=l+2)
      mc[99] = ((-3+2*ll-4*mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[99] = 0;
    if (abs(m-1)<=l+4)
      mc[100] = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[100] = 0;
    if (abs(m-3)<=l+4)
      mc[101] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[101] = 0;
    if (abs(m+1)<=l+4)
      mc[102] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[102] = 0;
    if (abs(m+3)<=l+4)
      mc[103] = -(sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[103] = 0;
    if (abs(m-4)<=l-2)
      mc[104] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-3+2*ll)*sqrt(1+2*ll)*(60-104*ll-48*pow(ll,2)+32*pow(ll,3)));
    else
      mc[104] = 0;
    if (abs(m)<=l-2)
      mc[105] = -(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-4-ll+pow(ll,2)+pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[105] = 0;
    if (abs(m+4)<=l-2)
      mc[106] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(-3+2*ll)*sqrt(1+2*ll)*(60-104*ll-48*pow(ll,2)+32*pow(ll,3)));
    else
      mc[106] = 0;
    if (abs(m-4)<=l-4)
      mc[107] = (sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(240-736*ll+576*pow(ll,2)-128*pow(ll,3)));
    else
      mc[107] = 0;
    if (abs(m)<=l-4)
      mc[108] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[108] = 0;
    if (abs(m+4)<=l-4)
      mc[109] = (sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(240-736*ll+576*pow(ll,2)-128*pow(ll,3)));
    else
      mc[109] = 0;
    if (abs(m-4)<=l)
      mc[110] = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[110] = 0;
    if (abs(m)<=l)
      mc[111] = (6*pow(ll,3)+3*pow(ll,4)+2*ll*(-7+pow(mm,2))+pow(ll,2)*(-11+2*pow(mm,2))+3*(4-5*pow(mm,2)+pow(mm,4)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[111] = 0;
    if (abs(m+4)<=l)
      mc[112] = (-3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[112] = 0;
    if (abs(m-4)<=l+2)
      mc[113] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[113] = 0;
    if (abs(m)<=l+2)
      mc[114] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-2+3*ll+pow(ll,2)+pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[114] = 0;
    if (abs(m+4)<=l+2)
      mc[115] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[115] = 0;
    if (abs(m-4)<=l+4)
      mc[116] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[116] = 0;
    if (abs(m)<=l+4)
      mc[117] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[117] = 0;
    if (abs(m+4)<=l+4)
      mc[118] = -(sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[118] = 0;
    if (abs(m-1)<=l-2)
      mc[119] = (sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)-3*mm+4*pow(mm,2)-ll*(1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[119] = 0;
    if (abs(m-3)<=l-2)
      mc[120] = -((5+2*ll-4*mm)*sqrt(1+ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[120] = 0;
    if (abs(m+1)<=l-2)
      mc[121] = (sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)+3*mm+4*pow(mm,2)+ll*(-1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[121] = 0;
    if (abs(m+3)<=l-2)
      mc[122] = -(sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*(5+2*ll+4*mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[122] = 0;
    if (abs(m-1)<=l-4)
      mc[123] = (sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[123] = 0;
    if (abs(m-3)<=l-4)
      mc[124] = (sqrt(ll-mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[124] = 0;
    if (abs(m+1)<=l-4)
      mc[125] = (sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[125] = 0;
    if (abs(m+3)<=l-4)
      mc[126] = (sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[126] = 0;
    if (abs(m-1)<=l)
      mc[127] = (sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-6+ll+pow(ll,2)-3*mm+3*pow(mm,2)))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[127] = 0;
    if (abs(m-3)<=l)
      mc[128] = (3*(3-2*mm)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[128] = 0;
    if (abs(m+1)<=l)
      mc[129] = -(sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(ll+pow(ll,2)+3*(-2+mm+pow(mm,2))))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[129] = 0;
    if (abs(m+3)<=l)
      mc[130] = (3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[130] = 0;
    if (abs(m-1)<=l+2)
      mc[131] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(4-2*pow(ll,2)+mm-4*pow(mm,2)-ll*(5+2*mm)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[131] = 0;
    if (abs(m-3)<=l+2)
      mc[132] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(ll+mm)*(-3+2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[132] = 0;
    if (abs(m+1)<=l+2)
      mc[133] = -(sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(-4+2*pow(ll,2)+ll*(5-2*mm)+mm+4*pow(mm,2)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[133] = 0;
    if (abs(m+3)<=l+2)
      mc[134] = ((-3+2*ll-4*mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[134] = 0;
    if (abs(m-1)<=l+4)
      mc[135] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[135] = 0;
    if (abs(m-3)<=l+4)
      mc[136] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(1+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[136] = 0;
    if (abs(m+1)<=l+4)
      mc[137] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[137] = 0;
    if (abs(m+3)<=l+4)
      mc[138] = -(sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[138] = 0;
    if (abs(m-2)<=l-2)
      mc[139] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-5+4*ll*(-1+mm)+6*mm-4*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[139] = 0;
    if (abs(m)<=l-2)
      mc[140] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-1+4*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[140] = 0;
    if (abs(m+2)<=l-2)
      mc[141] = -(sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(5+6*mm+4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[141] = 0;
    if (abs(m-2)<=l-4)
      mc[142] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[142] = 0;
    if (abs(m)<=l-4)
      mc[143] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3)));
    else
      mc[143] = 0;
    if (abs(m+2)<=l-4)
      mc[144] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[144] = 0;
    if (abs(m-2)<=l)
      mc[145] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-6+ll+pow(ll,2)+6*mm-3*pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[145] = 0;
    if (abs(m)<=l)
      mc[146] = (3+2*pow(ll,3)+pow(ll,4)-3*pow(mm,4)+2*pow(ll,2)*(-2+pow(mm,2))+ll*(-5+2*pow(mm,2)))/(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4));
    else
      mc[146] = 0;
    if (abs(m+2)<=l)
      mc[147] = -(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(ll+pow(ll,2)-3*(2+2*mm+pow(mm,2))))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[147] = 0;
    if (abs(m-2)<=l+2)
      mc[148] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(1+4*ll*(-1+mm)-2*mm+4*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[148] = 0;
    if (abs(m)<=l+2)
      mc[149] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-1+4*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[149] = 0;
    if (abs(m+2)<=l+2)
      mc[150] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-1-2*mm-4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[150] = 0;
    if (abs(m-2)<=l+4)
      mc[151] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)));
    else
      mc[151] = 0;
    if (abs(m)<=l+4)
      mc[152] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3))));
    else
      mc[152] = 0;
    if (abs(m+2)<=l+4)
      mc[153] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)));
    else
      mc[153] = 0;
    if (abs(m-2)<=l-2)
      mc[154] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-2*pow(ll,2)+3*mm+2*ll*mm-2*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[154] = 0;
    if (abs(m-4)<=l-2)
      mc[155] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3)));
    else
      mc[155] = 0;
    if (abs(m+2)<=l-2)
      mc[156] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(-5+2*pow(ll,2)+3*mm+2*ll*mm+2*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[156] = 0;
    if (abs(m+4)<=l-2)
      mc[157] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(-3+2*ll)*sqrt(1+2*ll)*(60-104*ll-48*pow(ll,2)+32*pow(ll,3)));
    else
      mc[157] = 0;
    if (abs(m-2)<=l-4)
      mc[158] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[158] = 0;
    if (abs(m-4)<=l-4)
      mc[159] = (sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[159] = 0;
    if (abs(m+2)<=l-4)
      mc[160] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[160] = 0;
    if (abs(m+4)<=l-4)
      mc[161] = (sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(240-736*ll+576*pow(ll,2)-128*pow(ll,3)));
    else
      mc[161] = 0;
    if (abs(m-2)<=l)
      mc[162] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+ll+pow(ll,2)-2*mm+pow(mm,2)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[162] = 0;
    if (abs(m-4)<=l)
      mc[163] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[163] = 0;
    if (abs(m+2)<=l)
      mc[164] = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+ll+pow(ll,2)+2*mm+pow(mm,2)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[164] = 0;
    if (abs(m+4)<=l)
      mc[165] = (-3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[165] = 0;
    if (abs(m-2)<=l+2)
      mc[166] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(3-2*pow(ll,2)+mm-2*pow(mm,2)-2*ll*(2+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[166] = 0;
    if (abs(m-4)<=l+2)
      mc[167] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[167] = 0;
    if (abs(m+2)<=l+2)
      mc[168] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-3+2*pow(ll,2)-2*ll*(-2+mm)+mm+2*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[168] = 0;
    if (abs(m+4)<=l+2)
      mc[169] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[169] = 0;
    if (abs(m-2)<=l+4)
      mc[170] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[170] = 0;
    if (abs(m-4)<=l+4)
      mc[171] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[171] = 0;
    if (abs(m+2)<=l+4)
      mc[172] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[172] = 0;
    if (abs(m+4)<=l+4)
      mc[173] = -(sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[173] = 0;
    if (abs(m-1)<=l-2)
      mc[174] = (sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(7+ll-2*pow(ll,2)+3*mm+2*ll*mm-4*pow(mm,2)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[174] = 0;
    if (abs(m-3)<=l-2)
      mc[175] = -((5+2*ll-4*mm)*sqrt(1+ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[175] = 0;
    if (abs(m+1)<=l-2)
      mc[176] = (sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)+3*mm+4*pow(mm,2)+ll*(-1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[176] = 0;
    if (abs(m+3)<=l-2)
      mc[177] = (sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*(5+2*ll+4*mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[177] = 0;
    if (abs(m-1)<=l-4)
      mc[178] = (sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[178] = 0;
    if (abs(m-3)<=l-4)
      mc[179] = (sqrt(ll-mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[179] = 0;
    if (abs(m+1)<=l-4)
      mc[180] = (sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[180] = 0;
    if (abs(m+3)<=l-4)
      mc[181] = (sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[181] = 0;
    if (abs(m-1)<=l)
      mc[182] = -(sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-6+ll+pow(ll,2)-3*mm+3*pow(mm,2)))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[182] = 0;
    if (abs(m-3)<=l)
      mc[183] = (3*(3-2*mm)*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[183] = 0;
    if (abs(m+1)<=l)
      mc[184] = -(sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(ll+pow(ll,2)+3*(-2+mm+pow(mm,2))))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[184] = 0;
    if (abs(m+3)<=l)
      mc[185] = (-3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[185] = 0;
    if (abs(m-1)<=l+2)
      mc[186] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-4+2*pow(ll,2)-mm+4*pow(mm,2)+ll*(5+2*mm)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[186] = 0;
    if (abs(m-3)<=l+2)
      mc[187] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(ll+mm)*(-3+2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[187] = 0;
    if (abs(m+1)<=l+2)
      mc[188] = -(sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(-4+2*pow(ll,2)+ll*(5-2*mm)+mm+4*pow(mm,2)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[188] = 0;
    if (abs(m+3)<=l+2)
      mc[189] = (sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*(3-2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[189] = 0;
    if (abs(m-1)<=l+4)
      mc[190] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[190] = 0;
    if (abs(m-3)<=l+4)
      mc[191] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(1+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[191] = 0;
    if (abs(m+1)<=l+4)
      mc[192] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[192] = 0;
    if (abs(m+3)<=l+4)
      mc[193] = (sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[193] = 0;
    if (abs(m-2)<=l-2)
      mc[194] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-4*ll*(-1+mm)-6*mm+4*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[194] = 0;
    if (abs(m+2)<=l-2)
      mc[195] = -(sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(5+6*mm+4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[195] = 0;
    if (abs(m-2)<=l-4)
      mc[196] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3)));
    else
      mc[196] = 0;
    if (abs(m+2)<=l-4)
      mc[197] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[197] = 0;
    if (abs(m-2)<=l)
      mc[198] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-6+ll+pow(ll,2)+6*mm-3*pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[198] = 0;
    if (abs(m+2)<=l)
      mc[199] = -(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(ll+pow(ll,2)-3*(2+2*mm+pow(mm,2))))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[199] = 0;
    if (abs(m-2)<=l+2)
      mc[200] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(1+4*ll*(-1+mm)-2*mm+4*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[200] = 0;
    if (abs(m+2)<=l+2)
      mc[201] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-1-2*mm-4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[201] = 0;
    if (abs(m-2)<=l+4)
      mc[202] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))));
    else
      mc[202] = 0;
    if (abs(m+2)<=l+4)
      mc[203] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)));
    else
      mc[203] = 0;
    if (abs(m-1)<=l-2)
      mc[204] = (sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(8-2*pow(ll,2)+ll*(3-2*mm)-3*mm+4*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[204] = 0;
    if (abs(m+1)<=l-2)
      mc[205] = (sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-8+2*pow(ll,2)-3*mm-4*pow(mm,2)-ll*(3+2*mm)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[205] = 0;
    if (abs(m-1)<=l-4)
      mc[206] = (sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3)));
    else
      mc[206] = 0;
    if (abs(m+1)<=l-4)
      mc[207] = -((sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3))));
    else
      mc[207] = 0;
    if (abs(m-1)<=l)
      mc[208] = (-3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-3+ll+pow(ll,2)+mm-pow(mm,2)))/(2.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[208] = 0;
    if (abs(m+1)<=l)
      mc[209] = (-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(-3+ll+pow(ll,2)-mm-pow(mm,2)))/(2.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[209] = 0;
    if (abs(m-1)<=l+2)
      mc[210] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-3+2*pow(ll,2)+ll*(7-2*mm)+mm-4*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[210] = 0;
    if (abs(m+1)<=l+2)
      mc[211] = (sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3-2*pow(ll,2)+mm+4*pow(mm,2)-ll*(7+2*mm)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[211] = 0;
    if (abs(m-1)<=l+4)
      mc[212] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3)));
    else
      mc[212] = 0;
    if (abs(m+1)<=l+4)
      mc[213] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3))));
    else
      mc[213] = 0;
    if (abs(m-2)<=l-2)
      mc[214] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-2*pow(ll,2)+3*mm+2*ll*mm-2*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[214] = 0;
    if (abs(m-4)<=l-2)
      mc[215] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3)));
    else
      mc[215] = 0;
    if (abs(m)<=l-2)
      mc[216] = (-3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-4-ll+pow(ll,2)+pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[216] = 0;
    if (abs(m+2)<=l-2)
      mc[217] = -(sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(-5+2*pow(ll,2)+3*mm+2*ll*mm+2*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[217] = 0;
    if (abs(m+4)<=l-2)
      mc[218] = -(sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(4.*sqrt(-3+2*ll)*sqrt(1+2*ll)*(15-26*ll-12*pow(ll,2)+8*pow(ll,3)));
    else
      mc[218] = 0;
    if (abs(m-2)<=l-4)
      mc[219] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[219] = 0;
    if (abs(m-4)<=l-4)
      mc[220] = (sqrt(-7+ll+mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[220] = 0;
    if (abs(m)<=l-4)
      mc[221] = (3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[221] = 0;
    if (abs(m+2)<=l-4)
      mc[222] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[222] = 0;
    if (abs(m+4)<=l-4)
      mc[223] = (sqrt(-7+ll-mm)*sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm))/(16.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[223] = 0;
    if (abs(m-2)<=l)
      mc[224] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+ll+pow(ll,2)-2*mm+pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[224] = 0;
    if (abs(m-4)<=l)
      mc[225] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[225] = 0;
    if (abs(m)<=l)
      mc[226] = (18*pow(ll,3)+9*pow(ll,4)+6*ll*(-7+pow(mm,2))+pow(ll,2)*(-33+6*pow(mm,2))+9*(4-5*pow(mm,2)+pow(mm,4)))/(4.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[226] = 0;
    if (abs(m+2)<=l)
      mc[227] = (3*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+ll+pow(ll,2)+2*mm+pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[227] = 0;
    if (abs(m+4)<=l)
      mc[228] = (3*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[228] = 0;
    if (abs(m-2)<=l+2)
      mc[229] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(3-2*pow(ll,2)+mm-2*pow(mm,2)-2*ll*(2+mm)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[229] = 0;
    if (abs(m-4)<=l+2)
      mc[230] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[230] = 0;
    if (abs(m)<=l+2)
      mc[231] = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-2+3*ll+pow(ll,2)+pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[231] = 0;
    if (abs(m+2)<=l+2)
      mc[232] = -(sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(-3+2*pow(ll,2)-2*ll*(-2+mm)+mm+2*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[232] = 0;
    if (abs(m+4)<=l+2)
      mc[233] = -(sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(4.*sqrt(1+2*ll)*sqrt(5+2*ll)*(-21+22*ll+36*pow(ll,2)+8*pow(ll,3)));
    else
      mc[233] = 0;
    if (abs(m-2)<=l+4)
      mc[234] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)));
    else
      mc[234] = 0;
    if (abs(m-4)<=l+4)
      mc[235] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(8+ll-mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[235] = 0;
    if (abs(m)<=l+4)
      mc[236] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[236] = 0;
    if (abs(m+2)<=l+4)
      mc[237] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3)));
    else
      mc[237] = 0;
    if (abs(m+4)<=l+4)
      mc[238] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm)*sqrt(8+ll+mm))/(16.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[238] = 0;
    if (abs(m-1)<=l-2)
      mc[239] = (3*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)-3*mm+4*pow(mm,2)-ll*(1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[239] = 0;
    if (abs(m-3)<=l-2)
      mc[240] = ((5+2*ll-4*mm)*sqrt(1+ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[240] = 0;
    if (abs(m+1)<=l-2)
      mc[241] = (3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-7+2*pow(ll,2)+3*mm+4*pow(mm,2)+ll*(-1+2*mm)))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[241] = 0;
    if (abs(m+3)<=l-2)
      mc[242] = (sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*(5+2*ll+4*mm))/(8.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[242] = 0;
    if (abs(m-1)<=l-4)
      mc[243] = (-3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[243] = 0;
    if (abs(m-3)<=l-4)
      mc[244] = (sqrt(ll-mm)*sqrt(-6+ll+mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[244] = 0;
    if (abs(m+1)<=l-4)
      mc[245] = (-3*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(8.*sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[245] = 0;
    if (abs(m+3)<=l-4)
      mc[246] = (sqrt(-6+ll-mm)*sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(120-368*ll+288*pow(ll,2)-64*pow(ll,3)));
    else
      mc[246] = 0;
    if (abs(m-1)<=l)
      mc[247] = (3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-6+ll+pow(ll,2)-3*mm+3*pow(mm,2)))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[247] = 0;
    if (abs(m-3)<=l)
      mc[248] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[248] = 0;
    if (abs(m+1)<=l)
      mc[249] = (-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(ll+pow(ll,2)+3*(-2+mm+pow(mm,2))))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[249] = 0;
    if (abs(m+3)<=l)
      mc[250] = (-3*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3+2*mm))/(8.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[250] = 0;
    if (abs(m-1)<=l+2)
      mc[251] = (-3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-4+2*pow(ll,2)-mm+4*pow(mm,2)+ll*(5+2*mm)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[251] = 0;
    if (abs(m-3)<=l+2)
      mc[252] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(ll+mm)*(-3+2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[252] = 0;
    if (abs(m+1)<=l+2)
      mc[253] = (-3*sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(-4+2*pow(ll,2)+ll*(5-2*mm)+mm+4*pow(mm,2)))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[253] = 0;
    if (abs(m+3)<=l+2)
      mc[254] = (sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*(3-2*ll+4*mm))/(8.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[254] = 0;
    if (abs(m-1)<=l+4)
      mc[255] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[255] = 0;
    if (abs(m-3)<=l+4)
      mc[256] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(7+ll-mm)*sqrt(1+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[256] = 0;
    if (abs(m+1)<=l+4)
      mc[257] = (3*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(8.*sqrt(1+2*ll)*(7+2*ll)*sqrt(9+2*ll)*(15+16*ll+4*pow(ll,2)));
    else
      mc[257] = 0;
    if (abs(m+3)<=l+4)
      mc[258] = (sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm)*sqrt(7+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(840+1136*ll+480*pow(ll,2)+64*pow(ll,3)));
    else
      mc[258] = 0;
    if (abs(m-2)<=l-2)
      mc[259] = (sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(5-4*ll*(-1+mm)-6*mm+4*pow(mm,2)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[259] = 0;
    if (abs(m)<=l-2)
      mc[260] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-1+4*pow(mm,2)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[260] = 0;
    if (abs(m+2)<=l-2)
      mc[261] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*(5+6*mm+4*pow(mm,2)+4*ll*(1+mm)))/(4.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[261] = 0;
    if (abs(m-2)<=l-4)
      mc[262] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-5+ll+mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3)));
    else
      mc[262] = 0;
    if (abs(m)<=l-4)
      mc[263] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3)));
    else
      mc[263] = 0;
    if (abs(m+2)<=l-4)
      mc[264] = (sqrt(-5+ll-mm)*sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(4.*sqrt(-7+2*ll)*sqrt(1+2*ll)*(15-46*ll+36*pow(ll,2)-8*pow(ll,3)));
    else
      mc[264] = 0;
    if (abs(m-2)<=l)
      mc[265] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-6+ll+pow(ll,2)+6*mm-3*pow(mm,2)))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[265] = 0;
    if (abs(m)<=l)
      mc[266] = (3+2*pow(ll,3)+pow(ll,4)-3*pow(mm,4)+2*pow(ll,2)*(-2+pow(mm,2))+ll*(-5+2*pow(mm,2)))/(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4));
    else
      mc[266] = 0;
    if (abs(m+2)<=l)
      mc[267] = (sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(ll+pow(ll,2)-3*(2+2*mm+pow(mm,2))))/(2.*(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4)));
    else
      mc[267] = 0;
    if (abs(m-2)<=l+2)
      mc[268] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*(1+4*ll*(-1+mm)-2*mm+4*pow(mm,2)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[268] = 0;
    if (abs(m)<=l+2)
      mc[269] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-1+4*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[269] = 0;
    if (abs(m+2)<=l+2)
      mc[270] = (sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*(1+2*mm+4*pow(mm,2)-4*ll*(1+mm)))/(4.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[270] = 0;
    if (abs(m-2)<=l+4)
      mc[271] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(6+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))));
    else
      mc[271] = 0;
    if (abs(m)<=l+4)
      mc[272] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3))));
    else
      mc[272] = 0;
    if (abs(m+2)<=l+4)
      mc[273] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm)*sqrt(6+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(420+568*ll+240*pow(ll,2)+32*pow(ll,3))));
    else
      mc[273] = 0;
    if (abs(m-1)<=l-2)
      mc[274] = (sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-8+2*pow(ll,2)+3*mm-4*pow(mm,2)+ll*(-3+2*mm)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[274] = 0;
    if (abs(m+1)<=l-2)
      mc[275] = (sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(ll+mm)*(-8+2*pow(ll,2)-3*mm-4*pow(mm,2)-ll*(3+2*mm)))/(2.*(-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[275] = 0;
    if (abs(m-1)<=l-4)
      mc[276] = -((sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-4+ll+mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3))));
    else
      mc[276] = 0;
    if (abs(m+1)<=l-4)
      mc[277] = -((sqrt(-4+ll-mm)*sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*sqrt(1+2*ll)*(30-92*ll+72*pow(ll,2)-16*pow(ll,3))));
    else
      mc[277] = 0;
    if (abs(m-1)<=l)
      mc[278] = (3*sqrt(1+ll-mm)*sqrt(ll+mm)*(-1+2*mm)*(-3+ll+pow(ll,2)+mm-pow(mm,2)))/(2.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[278] = 0;
    if (abs(m+1)<=l)
      mc[279] = (-3*sqrt(ll-mm)*sqrt(1+ll+mm)*(1+2*mm)*(-3+ll+pow(ll,2)-mm-pow(mm,2)))/(2.*(-15+4*ll+4*pow(ll,2))*(-3+4*ll+4*pow(ll,2)));
    else
      mc[279] = 0;
    if (abs(m-1)<=l+2)
      mc[280] = -(sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*(-3+2*pow(ll,2)+ll*(7-2*mm)+mm-4*pow(mm,2)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[280] = 0;
    if (abs(m+1)<=l+2)
      mc[281] = (sqrt(1+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*(3-2*pow(ll,2)+mm+4*pow(mm,2)-ll*(7+2*mm)))/(2.*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[281] = 0;
    if (abs(m-1)<=l+4)
      mc[282] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(5+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3))));
    else
      mc[282] = 0;
    if (abs(m+1)<=l+4)
      mc[283] = -((sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm)*sqrt(5+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(210+284*ll+120*pow(ll,2)+16*pow(ll,3))));
    else
      mc[283] = 0;
    if (abs(m)<=l-2)
      mc[284] = (2*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-1+ll+mm)*sqrt(ll+mm)*(-7-2*ll+2*pow(ll,2)-2*pow(mm,2)))/((-5+2*ll)*sqrt(-3+2*ll)*(-1+2*ll)*sqrt(1+2*ll)*(3+2*ll));
    else
      mc[284] = 0;
    if (abs(m)<=l-4)
      mc[285] = (sqrt(-3+ll-mm)*sqrt(-2+ll-mm)*sqrt(-1+ll-mm)*sqrt(ll-mm)*sqrt(-3+ll+mm)*sqrt(-2+ll+mm)*sqrt(-1+ll+mm)*sqrt(ll+mm))/(sqrt(-7+2*ll)*(-5+2*ll)*sqrt(1+2*ll)*(3-8*ll+4*pow(ll,2)));
    else
      mc[285] = 0;
    if (abs(m)<=l)
      mc[286] = (3*(3+4*pow(ll,3)+2*pow(ll,4)+10*pow(mm,2)+2*pow(mm,4)-4*ll*(2+pow(mm,2))-2*pow(ll,2)*(3+2*pow(mm,2))))/(45-72*ll-56*pow(ll,2)+32*pow(ll,3)+16*pow(ll,4));
    else
      mc[286] = 0;
    if (abs(m)<=l+2)
      mc[287] = (2*sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*(-3+6*ll+2*pow(ll,2)-2*pow(mm,2)))/((-1+2*ll)*sqrt(1+2*ll)*(3+2*ll)*sqrt(5+2*ll)*(7+2*ll));
    else
      mc[287] = 0;
    if (abs(m)<=l+4)
      mc[288] = (sqrt(1+ll-mm)*sqrt(2+ll-mm)*sqrt(3+ll-mm)*sqrt(4+ll-mm)*sqrt(1+ll+mm)*sqrt(2+ll+mm)*sqrt(3+ll+mm)*sqrt(4+ll+mm))/(sqrt(1+2*ll)*sqrt(9+2*ll)*(105+142*ll+60*pow(ll,2)+8*pow(ll,3)));
    else
      mc[288] = 0;
  }
}

static thread_return_t thread_function_0(void *arg) {
  struct job_info_s* job_info = (struct job_info_s*)arg;
  while (1) {
    REAL* psidot;
    REAL* psi;
    param_list_t *param;
    REAL a2[3][3];
    REAL a4[3][3][3][3];
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
    memcpy(&(a2),&(thread_pass_args_0.a2),sizeof(a2));
    memcpy(&(a4),&(thread_pass_args_0.a4),sizeof(a4));
    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {
      REAL *mc = mult_constant_0 + mc_count_0*mc_ind(l,m);
      REAL xx, xy, xz, yy, yz, zz, xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz;
      REAL temp1;
      xx = (mc[0]*psi[ind(l-2,m-2,0)]+mc[1]*psi[ind(l-2,m,0)]+mc[2]*psi[ind(l-2,m+2,0)]+mc[3]*psi[ind(l,m-2,0)]+mc[4]*psi[ind(l,m,0)]+mc[5]*psi[ind(l,m+2,0)]+mc[6]*psi[ind(l+2,m-2,0)]+mc[7]*psi[ind(l+2,m,0)]+mc[8]*psi[ind(l+2,m+2,0)]);
      xy = ((-mc[9])*psi[ind(l-2,m-2,1)]+(-mc[10])*psi[ind(l-2,m+2,1)]+(-mc[11])*psi[ind(l,m-2,1)]+(-mc[12])*psi[ind(l,m+2,1)]+(-mc[13])*psi[ind(l+2,m-2,1)]+(-mc[14])*psi[ind(l+2,m+2,1)]);
      xz = (mc[15]*psi[ind(l-2,m-1,0)]+mc[16]*psi[ind(l-2,m+1,0)]+mc[17]*psi[ind(l,m-1,0)]+mc[18]*psi[ind(l,m+1,0)]+mc[19]*psi[ind(l+2,m-1,0)]+mc[20]*psi[ind(l+2,m+1,0)]);
      yy = (mc[21]*psi[ind(l-2,m-2,0)]+mc[22]*psi[ind(l-2,m,0)]+mc[23]*psi[ind(l-2,m+2,0)]+mc[24]*psi[ind(l,m-2,0)]+mc[25]*psi[ind(l,m,0)]+mc[26]*psi[ind(l,m+2,0)]+mc[27]*psi[ind(l+2,m-2,0)]+mc[28]*psi[ind(l+2,m,0)]+mc[29]*psi[ind(l+2,m+2,0)]);
      yz = ((-mc[30])*psi[ind(l-2,m-1,1)]+(-mc[31])*psi[ind(l-2,m+1,1)]+(-mc[32])*psi[ind(l,m-1,1)]+(-mc[33])*psi[ind(l,m+1,1)]+(-mc[34])*psi[ind(l+2,m-1,1)]+(-mc[35])*psi[ind(l+2,m+1,1)]);
      zz = (mc[36]*psi[ind(l-2,m,0)]+mc[37]*psi[ind(l,m,0)]+mc[38]*psi[ind(l+2,m,0)]);
      xxxx = (mc[39]*psi[ind(l-2,m-2,0)]+mc[40]*psi[ind(l-2,m-4,0)]+mc[41]*psi[ind(l-2,m,0)]+mc[42]*psi[ind(l-2,m+2,0)]+mc[43]*psi[ind(l-2,m+4,0)]+mc[44]*psi[ind(l-4,m-2,0)]+mc[45]*psi[ind(l-4,m-4,0)]+mc[46]*psi[ind(l-4,m,0)]+mc[47]*psi[ind(l-4,m+2,0)]+mc[48]*psi[ind(l-4,m+4,0)]+mc[49]*psi[ind(l,m-2,0)]+mc[50]*psi[ind(l,m-4,0)]+mc[51]*psi[ind(l,m,0)]+mc[52]*psi[ind(l,m+2,0)]+mc[53]*psi[ind(l,m+4,0)]+mc[54]*psi[ind(l+2,m-2,0)]+mc[55]*psi[ind(l+2,m-4,0)]+mc[56]*psi[ind(l+2,m,0)]+mc[57]*psi[ind(l+2,m+2,0)]+mc[58]*psi[ind(l+2,m+4,0)]+mc[59]*psi[ind(l+4,m-2,0)]+mc[60]*psi[ind(l+4,m-4,0)]+mc[61]*psi[ind(l+4,m,0)]+mc[62]*psi[ind(l+4,m+2,0)]+mc[63]*psi[ind(l+4,m+4,0)]);
      xxxy = ((-mc[64])*psi[ind(l-2,m-2,1)]+(-mc[65])*psi[ind(l-2,m-4,1)]+(-mc[66])*psi[ind(l-2,m+2,1)]+(-mc[67])*psi[ind(l-2,m+4,1)]+(-mc[68])*psi[ind(l-4,m-2,1)]+(-mc[69])*psi[ind(l-4,m-4,1)]+(-mc[70])*psi[ind(l-4,m+2,1)]+(-mc[71])*psi[ind(l-4,m+4,1)]+(-mc[72])*psi[ind(l,m-2,1)]+(-mc[73])*psi[ind(l,m-4,1)]+(-mc[74])*psi[ind(l,m+2,1)]+(-mc[75])*psi[ind(l,m+4,1)]+(-mc[76])*psi[ind(l+2,m-2,1)]+(-mc[77])*psi[ind(l+2,m-4,1)]+(-mc[78])*psi[ind(l+2,m+2,1)]+(-mc[79])*psi[ind(l+2,m+4,1)]+(-mc[80])*psi[ind(l+4,m-2,1)]+(-mc[81])*psi[ind(l+4,m-4,1)]+(-mc[82])*psi[ind(l+4,m+2,1)]+(-mc[83])*psi[ind(l+4,m+4,1)]);
      xxxz = (mc[84]*psi[ind(l-2,m-1,0)]+mc[85]*psi[ind(l-2,m-3,0)]+mc[86]*psi[ind(l-2,m+1,0)]+mc[87]*psi[ind(l-2,m+3,0)]+mc[88]*psi[ind(l-4,m-1,0)]+mc[89]*psi[ind(l-4,m-3,0)]+mc[90]*psi[ind(l-4,m+1,0)]+mc[91]*psi[ind(l-4,m+3,0)]+mc[92]*psi[ind(l,m-1,0)]+mc[93]*psi[ind(l,m-3,0)]+mc[94]*psi[ind(l,m+1,0)]+mc[95]*psi[ind(l,m+3,0)]+mc[96]*psi[ind(l+2,m-1,0)]+mc[97]*psi[ind(l+2,m-3,0)]+mc[98]*psi[ind(l+2,m+1,0)]+mc[99]*psi[ind(l+2,m+3,0)]+mc[100]*psi[ind(l+4,m-1,0)]+mc[101]*psi[ind(l+4,m-3,0)]+mc[102]*psi[ind(l+4,m+1,0)]+mc[103]*psi[ind(l+4,m+3,0)]);
      xxyy = (mc[104]*psi[ind(l-2,m-4,0)]+mc[105]*psi[ind(l-2,m,0)]+mc[106]*psi[ind(l-2,m+4,0)]+mc[107]*psi[ind(l-4,m-4,0)]+mc[108]*psi[ind(l-4,m,0)]+mc[109]*psi[ind(l-4,m+4,0)]+mc[110]*psi[ind(l,m-4,0)]+mc[111]*psi[ind(l,m,0)]+mc[112]*psi[ind(l,m+4,0)]+mc[113]*psi[ind(l+2,m-4,0)]+mc[114]*psi[ind(l+2,m,0)]+mc[115]*psi[ind(l+2,m+4,0)]+mc[116]*psi[ind(l+4,m-4,0)]+mc[117]*psi[ind(l+4,m,0)]+mc[118]*psi[ind(l+4,m+4,0)]);
      xxyz = ((-mc[119])*psi[ind(l-2,m-1,1)]+(-mc[120])*psi[ind(l-2,m-3,1)]+(-mc[121])*psi[ind(l-2,m+1,1)]+(-mc[122])*psi[ind(l-2,m+3,1)]+(-mc[123])*psi[ind(l-4,m-1,1)]+(-mc[124])*psi[ind(l-4,m-3,1)]+(-mc[125])*psi[ind(l-4,m+1,1)]+(-mc[126])*psi[ind(l-4,m+3,1)]+(-mc[127])*psi[ind(l,m-1,1)]+(-mc[128])*psi[ind(l,m-3,1)]+(-mc[129])*psi[ind(l,m+1,1)]+(-mc[130])*psi[ind(l,m+3,1)]+(-mc[131])*psi[ind(l+2,m-1,1)]+(-mc[132])*psi[ind(l+2,m-3,1)]+(-mc[133])*psi[ind(l+2,m+1,1)]+(-mc[134])*psi[ind(l+2,m+3,1)]+(-mc[135])*psi[ind(l+4,m-1,1)]+(-mc[136])*psi[ind(l+4,m-3,1)]+(-mc[137])*psi[ind(l+4,m+1,1)]+(-mc[138])*psi[ind(l+4,m+3,1)]);
      xxzz = (mc[139]*psi[ind(l-2,m-2,0)]+mc[140]*psi[ind(l-2,m,0)]+mc[141]*psi[ind(l-2,m+2,0)]+mc[142]*psi[ind(l-4,m-2,0)]+mc[143]*psi[ind(l-4,m,0)]+mc[144]*psi[ind(l-4,m+2,0)]+mc[145]*psi[ind(l,m-2,0)]+mc[146]*psi[ind(l,m,0)]+mc[147]*psi[ind(l,m+2,0)]+mc[148]*psi[ind(l+2,m-2,0)]+mc[149]*psi[ind(l+2,m,0)]+mc[150]*psi[ind(l+2,m+2,0)]+mc[151]*psi[ind(l+4,m-2,0)]+mc[152]*psi[ind(l+4,m,0)]+mc[153]*psi[ind(l+4,m+2,0)]);
      xyyy = ((-mc[154])*psi[ind(l-2,m-2,1)]+(-mc[155])*psi[ind(l-2,m-4,1)]+(-mc[156])*psi[ind(l-2,m+2,1)]+(-mc[157])*psi[ind(l-2,m+4,1)]+(-mc[158])*psi[ind(l-4,m-2,1)]+(-mc[159])*psi[ind(l-4,m-4,1)]+(-mc[160])*psi[ind(l-4,m+2,1)]+(-mc[161])*psi[ind(l-4,m+4,1)]+(-mc[162])*psi[ind(l,m-2,1)]+(-mc[163])*psi[ind(l,m-4,1)]+(-mc[164])*psi[ind(l,m+2,1)]+(-mc[165])*psi[ind(l,m+4,1)]+(-mc[166])*psi[ind(l+2,m-2,1)]+(-mc[167])*psi[ind(l+2,m-4,1)]+(-mc[168])*psi[ind(l+2,m+2,1)]+(-mc[169])*psi[ind(l+2,m+4,1)]+(-mc[170])*psi[ind(l+4,m-2,1)]+(-mc[171])*psi[ind(l+4,m-4,1)]+(-mc[172])*psi[ind(l+4,m+2,1)]+(-mc[173])*psi[ind(l+4,m+4,1)]);
      xyyz = (mc[174]*psi[ind(l-2,m-1,0)]+mc[175]*psi[ind(l-2,m-3,0)]+mc[176]*psi[ind(l-2,m+1,0)]+mc[177]*psi[ind(l-2,m+3,0)]+mc[178]*psi[ind(l-4,m-1,0)]+mc[179]*psi[ind(l-4,m-3,0)]+mc[180]*psi[ind(l-4,m+1,0)]+mc[181]*psi[ind(l-4,m+3,0)]+mc[182]*psi[ind(l,m-1,0)]+mc[183]*psi[ind(l,m-3,0)]+mc[184]*psi[ind(l,m+1,0)]+mc[185]*psi[ind(l,m+3,0)]+mc[186]*psi[ind(l+2,m-1,0)]+mc[187]*psi[ind(l+2,m-3,0)]+mc[188]*psi[ind(l+2,m+1,0)]+mc[189]*psi[ind(l+2,m+3,0)]+mc[190]*psi[ind(l+4,m-1,0)]+mc[191]*psi[ind(l+4,m-3,0)]+mc[192]*psi[ind(l+4,m+1,0)]+mc[193]*psi[ind(l+4,m+3,0)]);
      xyzz = ((-mc[194])*psi[ind(l-2,m-2,1)]+(-mc[195])*psi[ind(l-2,m+2,1)]+(-mc[196])*psi[ind(l-4,m-2,1)]+(-mc[197])*psi[ind(l-4,m+2,1)]+(-mc[198])*psi[ind(l,m-2,1)]+(-mc[199])*psi[ind(l,m+2,1)]+(-mc[200])*psi[ind(l+2,m-2,1)]+(-mc[201])*psi[ind(l+2,m+2,1)]+(-mc[202])*psi[ind(l+4,m-2,1)]+(-mc[203])*psi[ind(l+4,m+2,1)]);
      xzzz = (mc[204]*psi[ind(l-2,m-1,0)]+mc[205]*psi[ind(l-2,m+1,0)]+mc[206]*psi[ind(l-4,m-1,0)]+mc[207]*psi[ind(l-4,m+1,0)]+mc[208]*psi[ind(l,m-1,0)]+mc[209]*psi[ind(l,m+1,0)]+mc[210]*psi[ind(l+2,m-1,0)]+mc[211]*psi[ind(l+2,m+1,0)]+mc[212]*psi[ind(l+4,m-1,0)]+mc[213]*psi[ind(l+4,m+1,0)]);
      yyyy = (mc[214]*psi[ind(l-2,m-2,0)]+mc[215]*psi[ind(l-2,m-4,0)]+mc[216]*psi[ind(l-2,m,0)]+mc[217]*psi[ind(l-2,m+2,0)]+mc[218]*psi[ind(l-2,m+4,0)]+mc[219]*psi[ind(l-4,m-2,0)]+mc[220]*psi[ind(l-4,m-4,0)]+mc[221]*psi[ind(l-4,m,0)]+mc[222]*psi[ind(l-4,m+2,0)]+mc[223]*psi[ind(l-4,m+4,0)]+mc[224]*psi[ind(l,m-2,0)]+mc[225]*psi[ind(l,m-4,0)]+mc[226]*psi[ind(l,m,0)]+mc[227]*psi[ind(l,m+2,0)]+mc[228]*psi[ind(l,m+4,0)]+mc[229]*psi[ind(l+2,m-2,0)]+mc[230]*psi[ind(l+2,m-4,0)]+mc[231]*psi[ind(l+2,m,0)]+mc[232]*psi[ind(l+2,m+2,0)]+mc[233]*psi[ind(l+2,m+4,0)]+mc[234]*psi[ind(l+4,m-2,0)]+mc[235]*psi[ind(l+4,m-4,0)]+mc[236]*psi[ind(l+4,m,0)]+mc[237]*psi[ind(l+4,m+2,0)]+mc[238]*psi[ind(l+4,m+4,0)]);
      yyyz = ((-mc[239])*psi[ind(l-2,m-1,1)]+(-mc[240])*psi[ind(l-2,m-3,1)]+(-mc[241])*psi[ind(l-2,m+1,1)]+(-mc[242])*psi[ind(l-2,m+3,1)]+(-mc[243])*psi[ind(l-4,m-1,1)]+(-mc[244])*psi[ind(l-4,m-3,1)]+(-mc[245])*psi[ind(l-4,m+1,1)]+(-mc[246])*psi[ind(l-4,m+3,1)]+(-mc[247])*psi[ind(l,m-1,1)]+(-mc[248])*psi[ind(l,m-3,1)]+(-mc[249])*psi[ind(l,m+1,1)]+(-mc[250])*psi[ind(l,m+3,1)]+(-mc[251])*psi[ind(l+2,m-1,1)]+(-mc[252])*psi[ind(l+2,m-3,1)]+(-mc[253])*psi[ind(l+2,m+1,1)]+(-mc[254])*psi[ind(l+2,m+3,1)]+(-mc[255])*psi[ind(l+4,m-1,1)]+(-mc[256])*psi[ind(l+4,m-3,1)]+(-mc[257])*psi[ind(l+4,m+1,1)]+(-mc[258])*psi[ind(l+4,m+3,1)]);
      yyzz = (mc[259]*psi[ind(l-2,m-2,0)]+mc[260]*psi[ind(l-2,m,0)]+mc[261]*psi[ind(l-2,m+2,0)]+mc[262]*psi[ind(l-4,m-2,0)]+mc[263]*psi[ind(l-4,m,0)]+mc[264]*psi[ind(l-4,m+2,0)]+mc[265]*psi[ind(l,m-2,0)]+mc[266]*psi[ind(l,m,0)]+mc[267]*psi[ind(l,m+2,0)]+mc[268]*psi[ind(l+2,m-2,0)]+mc[269]*psi[ind(l+2,m,0)]+mc[270]*psi[ind(l+2,m+2,0)]+mc[271]*psi[ind(l+4,m-2,0)]+mc[272]*psi[ind(l+4,m,0)]+mc[273]*psi[ind(l+4,m+2,0)]);
      yzzz = ((-mc[274])*psi[ind(l-2,m-1,1)]+(-mc[275])*psi[ind(l-2,m+1,1)]+(-mc[276])*psi[ind(l-4,m-1,1)]+(-mc[277])*psi[ind(l-4,m+1,1)]+(-mc[278])*psi[ind(l,m-1,1)]+(-mc[279])*psi[ind(l,m+1,1)]+(-mc[280])*psi[ind(l+2,m-1,1)]+(-mc[281])*psi[ind(l+2,m+1,1)]+(-mc[282])*psi[ind(l+4,m-1,1)]+(-mc[283])*psi[ind(l+4,m+1,1)]);
      zzzz = (mc[284]*psi[ind(l-2,m,0)]+mc[285]*psi[ind(l-4,m,0)]+mc[286]*psi[ind(l,m,0)]+mc[287]*psi[ind(l+2,m,0)]+mc[288]*psi[ind(l+4,m,0)]);
      temp1 = xxzz*pow(param->gamm[0*3+1],2)*a2[0][0] - 
              2*xxyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][0] + 
              xxyy*pow(param->gamm[0*3+2],2)*a2[0][0] + 
              2*xyzz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[0][0] - 
              2*xyyz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][0] + 
              yyzz*pow(param->gamm[1*3+1],2)*a2[0][0] - 
              2*xyyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][0] + 
              2*xzzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][0] + 
              2*xyyy*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][0] - 
              2*xyzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][0] - 
              2*yyyz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][0] + 
              2*yzzz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][0] + 
              yyyy*pow(param->gamm[1*3+2],2)*a2[0][0] - 
              2*yyzz*pow(param->gamm[1*3+2],2)*a2[0][0] + 
              zzzz*pow(param->gamm[1*3+2],2)*a2[0][0] - 
              2*xyzz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][0] + 
              2*xyyz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[0][0] - 
              2*yyzz*param->gamm[1*3+1]*param->gamm[2*3+2]*a2[0][0] + 
              2*yyyz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][0] - 
              2*yzzz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][0] + 
              yyzz*pow(param->gamm[2*3+2],2)*a2[0][0] - 
              2*xxzz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[0][1] - 
              2*xyzz*pow(param->gamm[0*3+1],2)*a2[0][1] + 
              2*xxyz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[0][1] + 
              2*xxxz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][1] + 
              2*xyyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][1] - 
              2*xzzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][1] - 
              2*xxxy*pow(param->gamm[0*3+2],2)*a2[0][1] + 
              2*xyzz*pow(param->gamm[0*3+2],2)*a2[0][1] - 
              2*xyzz*param->gamm[0*3+0]*param->gamm[1*3+1]*a2[0][1] - 
              2*yyzz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[0][1] + 
              2*xxyz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][1] - 
              2*yzzz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][1] + 
              2*xyyz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[0][1] - 
              2*xzzz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[0][1] + 
              2*xxyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][1] + 
              2*yyyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][1] - 
              2*yzzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][1] - 
              4*xxyy*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][1] + 
              2*xxzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][1] + 
              2*yyzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][1] - 
              2*zzzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][1] + 
              2*xyyz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][1] - 
              2*xyyy*pow(param->gamm[1*3+2],2)*a2[0][1] + 
              2*xyzz*pow(param->gamm[1*3+2],2)*a2[0][1] + 
              2*xyzz*param->gamm[0*3+0]*param->gamm[2*3+2]*a2[0][1] + 
              2*xxzz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][1] + 
              2*yyzz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][1] - 
              4*xxyz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[0][1] + 
              2*yzzz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[0][1] + 
              2*xyzz*param->gamm[1*3+1]*param->gamm[2*3+2]*a2[0][1] - 
              4*xyyz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][1] + 
              2*xzzz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][1] - 
              2*xyzz*pow(param->gamm[2*3+2],2)*a2[0][1] + 
              2*xxyz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[0][2] - 
              2*xxxz*pow(param->gamm[0*3+1],2)*a2[0][2] + 
              2*xyyz*pow(param->gamm[0*3+1],2)*a2[0][2] - 
              2*xxyy*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[0][2] + 
              2*xxxy*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][2] - 
              2*xyyy*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][2] + 
              2*xyzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][2] - 
              2*xyyz*pow(param->gamm[0*3+2],2)*a2[0][2] + 
              2*xyyz*param->gamm[0*3+0]*param->gamm[1*3+1]*a2[0][2] - 
              4*xxyz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[0][2] + 
              2*yyyz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[0][2] + 
              2*xxyy*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][2] + 
              2*yyzz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][2] - 
              2*xyyz*pow(param->gamm[1*3+1],2)*a2[0][2] - 
              2*xyyy*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[0][2] + 
              2*xyzz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[0][2] + 
              2*xxyy*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][2] - 
              4*xxzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][2] - 
              2*yyyy*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][2] + 
              2*yyzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][2] + 
              2*xxyz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][2] - 
              2*yyyz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][2] + 
              2*yzzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][2] + 
              2*xyyy*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][2] - 
              4*xyzz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][2] + 
              2*xyyz*pow(param->gamm[1*3+2],2)*a2[0][2] - 
              2*xzzz*pow(param->gamm[1*3+2],2)*a2[0][2] - 
              2*xyyz*param->gamm[0*3+0]*param->gamm[2*3+2]*a2[0][2] + 
              2*xxyz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][2] - 
              2*yyyz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][2] - 
              2*yyzz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[0][2] + 
              2*xyyz*param->gamm[1*3+1]*param->gamm[2*3+2]*a2[0][2] + 
              2*xyzz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][2] + 
              xxzz*pow(param->gamm[0*3+0],2)*a2[1][1] + 
              2*xyzz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[1][1] + 
              yyzz*pow(param->gamm[0*3+1],2)*a2[1][1] - 
              2*xxxz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[1][1] + 
              2*xzzz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[1][1] - 
              2*xxyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][1] + 
              2*yzzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][1] + 
              xxxx*pow(param->gamm[0*3+2],2)*a2[1][1] - 
              2*xxzz*pow(param->gamm[0*3+2],2)*a2[1][1] + 
              zzzz*pow(param->gamm[0*3+2],2)*a2[1][1] - 
              2*xxyz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[1][1] - 
              2*xyyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[1][1] + 
              2*xxxy*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][1] - 
              2*xyzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][1] + 
              xxyy*pow(param->gamm[1*3+2],2)*a2[1][1] - 
              2*xxzz*param->gamm[0*3+0]*param->gamm[2*3+2]*a2[1][1] - 
              2*xyzz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[1][1] + 
              2*xxxz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[1][1] - 
              2*xzzz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[1][1] + 
              2*xxyz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[1][1] + 
              xxzz*pow(param->gamm[2*3+2],2)*a2[1][1] - 
              2*xxyz*pow(param->gamm[0*3+0],2)*a2[1][2] + 
              2*xxxz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[1][2] - 
              4*xyyz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[1][2] + 
              2*xxyz*pow(param->gamm[0*3+1],2)*a2[1][2] - 
              2*yyyz*pow(param->gamm[0*3+1],2)*a2[1][2] + 
              2*xxxy*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[1][2] - 
              4*xyzz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[1][2] - 
              2*xxxx*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][2] + 
              2*xxyy*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][2] + 
              2*xxzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][2] - 
              4*yyzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][2] + 
              2*xxyz*pow(param->gamm[0*3+2],2)*a2[1][2] - 
              2*yzzz*pow(param->gamm[0*3+2],2)*a2[1][2] + 
              2*xxyz*param->gamm[0*3+0]*param->gamm[1*3+1]*a2[1][2] + 
              2*xyyz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[1][2] - 
              2*xxxy*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[1][2] + 
              2*xyzz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[1][2] + 
              2*xxyy*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[1][2] + 
              2*xxzz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[1][2] - 
              2*xxxy*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[1][2] + 
              2*xyyy*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[1][2] + 
              2*xyzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[1][2] - 
              2*xxxz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][2] + 
              2*xyyz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][2] + 
              2*xzzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][2] - 
              2*xxyy*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[1][2] - 
              2*xxyz*pow(param->gamm[1*3+2],2)*a2[1][2] + 
              2*xxyz*param->gamm[0*3+0]*param->gamm[2*3+2]*a2[1][2] - 
              2*xxxz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[1][2] + 
              2*xyyz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[1][2] + 
              2*xyzz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[1][2] - 
              2*xxyz*param->gamm[1*3+1]*param->gamm[2*3+2]*a2[1][2] - 
              2*xxzz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[1][2] + 
              xxyy*pow(param->gamm[0*3+0],2)*a2[2][2] - 
              2*xxxy*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[2][2] + 
              2*xyyy*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[2][2] + 
              xxxx*pow(param->gamm[0*3+1],2)*a2[2][2] - 
              2*xxyy*pow(param->gamm[0*3+1],2)*a2[2][2] + 
              yyyy*pow(param->gamm[0*3+1],2)*a2[2][2] + 
              2*xyyz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[2][2] - 
              2*xxyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[2][2] + 
              2*yyyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[2][2] + 
              yyzz*pow(param->gamm[0*3+2],2)*a2[2][2] - 
              2*xxyy*param->gamm[0*3+0]*param->gamm[1*3+1]*a2[2][2] + 
              2*xxxy*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[2][2] - 
              2*xyyy*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[2][2] - 
              2*xyyz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[2][2] + 
              xxyy*pow(param->gamm[1*3+1],2)*a2[2][2] - 
              2*xxyz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[2][2] + 
              2*xxxz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[2][2] - 
              2*xyyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[2][2] - 
              2*xyzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[2][2] + 
              2*xxyz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[2][2] + 
              xxzz*pow(param->gamm[1*3+2],2)*a2[2][2] + 
              zz*pow(param->gamm[0*3+1],2)*a4[0][0][0][0] - 
              2*yz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][0][0] + 
              yy*pow(param->gamm[0*3+2],2)*a4[0][0][0][0] - 
              2*zz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][0][0][1] + 
              2*yz*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][0][0][1] + 
              2*xz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][0][1] - 
              2*xy*pow(param->gamm[0*3+2],2)*a4[0][0][0][1] + 
              2*zz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][0][0][1] - 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][0][0][1] - 
              2*yz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][0][1] + 
              2*yy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][0][1] + 
              2*yz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][0][0][2] - 
              2*xz*pow(param->gamm[0*3+1],2)*a4[0][0][0][2] - 
              2*yy*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][0][0][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][0][2] + 
              2*zz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][0][2] - 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][0][2] - 
              2*yz*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][0][0][2] + 
              2*yy*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][0][0][2] + 
              zz*pow(param->gamm[0*3+0],2)*a4[0][0][1][1] - 
              2*zz*pow(param->gamm[0*3+1],2)*a4[0][0][1][1] - 
              2*xz*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][0][1][1] + 
              2*yz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][1][1] + 
              xx*pow(param->gamm[0*3+2],2)*a4[0][0][1][1] - 
              2*zz*param->gamm[0*3+0]*param->gamm[1*3+1]*a4[0][0][1][1] + 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][0][1][1] + 
              zz*pow(param->gamm[1*3+1],2)*a4[0][0][1][1] + 
              2*yz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][0][1][1] + 
              2*xz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][1][1] - 
              4*xy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][1][1] - 
              2*yz*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][0][1][1] + 
              yy*pow(param->gamm[1*3+2],2)*a4[0][0][1][1] - 
              2*yz*pow(param->gamm[0*3+0],2)*a4[0][0][1][2] + 
              2*xz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][0][1][2] + 
              2*yz*pow(param->gamm[0*3+1],2)*a4[0][0][1][2] + 
              2*xy*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][0][1][2] - 
              2*xx*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][1][2] - 
              2*yy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][1][2] - 
              2*zz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][1][2] + 
              2*yz*pow(param->gamm[0*3+2],2)*a4[0][0][1][2] + 
              2*yz*param->gamm[0*3+0]*param->gamm[1*3+1]*a4[0][0][1][2] - 
              4*xz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][0][1][2] + 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][0][1][2] - 
              2*yy*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][0][1][2] - 
              2*zz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][0][1][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][1][2] + 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][1][2] + 
              2*zz*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][0][1][2] - 
              2*yz*pow(param->gamm[1*3+2],2)*a4[0][0][1][2] + 
              2*yz*param->gamm[0*3+0]*param->gamm[2*3+2]*a4[0][0][1][2] + 
              2*xz*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][0][1][2] - 
              4*xy*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][0][1][2] - 
              2*yz*param->gamm[1*3+1]*param->gamm[2*3+2]*a4[0][0][1][2] + 
              2*yy*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][0][1][2] + 
              yy*pow(param->gamm[0*3+0],2)*a4[0][0][2][2] - 
              2*xy*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][0][2][2] + 
              xx*pow(param->gamm[0*3+1],2)*a4[0][0][2][2] + 
              2*yz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][2][2] - 
              2*yy*pow(param->gamm[0*3+2],2)*a4[0][0][2][2] + 
              2*yz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][0][2][2] - 
              4*xz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][2][2] + 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][2][2] + 
              zz*pow(param->gamm[1*3+2],2)*a4[0][0][2][2] - 
              2*yy*param->gamm[0*3+0]*param->gamm[2*3+2]*a4[0][0][2][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][0][2][2] - 
              2*yz*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][0][2][2] + 
              yy*pow(param->gamm[2*3+2],2)*a4[0][0][2][2] + 
              2*zz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][1][1][1] - 
              2*xz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][1][1][1] - 
              2*zz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][1][1][1] - 
              2*xz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][1][1][1] + 
              2*yz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][1][1] + 
              2*xx*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][1][1] + 
              2*xz*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][1][1][1] - 
              2*xy*pow(param->gamm[1*3+2],2)*a4[0][1][1][1] - 
              4*yz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][1][1][2] + 
              2*xz*pow(param->gamm[0*3+1],2)*a4[0][1][1][2] + 
              2*zz*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][1][1][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][1][1][2] - 
              2*xz*pow(param->gamm[0*3+2],2)*a4[0][1][1][2] + 
              2*xz*param->gamm[0*3+0]*param->gamm[1*3+1]*a4[0][1][1][2] + 
              2*yz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][1][1][2] - 
              2*xx*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][1][1][2] - 
              2*zz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][1][1][2] - 
              2*xz*pow(param->gamm[1*3+1],2)*a4[0][1][1][2] + 
              2*xy*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][1][1][2] - 
              2*xx*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][1][2] - 
              2*yy*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][1][2] - 
              2*zz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][1][2] + 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][1][2] + 
              2*xy*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][1][1][2] + 
              2*xz*pow(param->gamm[1*3+2],2)*a4[0][1][1][2] - 
              2*xz*param->gamm[0*3+0]*param->gamm[2*3+2]*a4[0][1][1][2] + 
              2*yz*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][1][1][2] + 
              2*xx*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][1][1][2] + 
              2*xz*param->gamm[1*3+1]*param->gamm[2*3+2]*a4[0][1][1][2] - 
              4*xy*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][1][1][2] + 
              2*yy*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][1][2][2] - 
              2*xy*pow(param->gamm[0*3+1],2)*a4[0][1][2][2] - 
              4*yz*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][1][2][2] + 
              2*xz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][1][2][2] + 
              2*xy*pow(param->gamm[0*3+2],2)*a4[0][1][2][2] - 
              2*xy*param->gamm[0*3+0]*param->gamm[1*3+1]*a4[0][1][2][2] + 
              2*xx*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][1][2][2] + 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][1][2][2] + 
              2*xz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][1][2][2] + 
              2*yz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][2][2] - 
              2*xx*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][2][2] - 
              2*yy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][2][2] - 
              2*zz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][2][2] - 
              4*xz*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][1][2][2] + 
              2*xy*pow(param->gamm[1*3+2],2)*a4[0][1][2][2] + 
              2*xy*param->gamm[0*3+0]*param->gamm[2*3+2]*a4[0][1][2][2] - 
              2*xx*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][1][2][2] - 
              2*yy*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][1][2][2] + 
              2*yz*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][1][2][2] + 
              2*xy*param->gamm[1*3+1]*param->gamm[2*3+2]*a4[0][1][2][2] + 
              2*xz*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][1][2][2] - 
              2*xy*pow(param->gamm[2*3+2],2)*a4[0][1][2][2] + 
              2*yy*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][2][2][2] - 
              2*xy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][2][2][2] - 
              2*xy*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][2][2][2] + 
              2*xx*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][2][2][2] + 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][2][2][2] - 
              2*xz*pow(param->gamm[1*3+2],2)*a4[0][2][2][2] - 
              2*yy*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][2][2][2] + 
              2*xy*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][2][2][2] + 
              zz*pow(param->gamm[0*3+1],2)*a4[1][1][1][1] - 
              2*xz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[1][1][1][1] + 
              xx*pow(param->gamm[1*3+2],2)*a4[1][1][1][1] - 
              2*yz*pow(param->gamm[0*3+1],2)*a4[1][1][1][2] + 
              2*zz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[1][1][1][2] + 
              2*xz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[1][1][1][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[1][1][1][2] - 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[1][1][1][2] - 
              2*xx*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[1][1][1][2] - 
              2*xz*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[1][1][1][2] + 
              2*xx*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[1][1][1][2] + 
              yy*pow(param->gamm[0*3+1],2)*a4[1][1][2][2] - 
              4*yz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[1][1][2][2] + 
              zz*pow(param->gamm[0*3+2],2)*a4[1][1][2][2] - 
              2*xy*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[1][1][2][2] + 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[1][1][2][2] + 
              xx*pow(param->gamm[1*3+1],2)*a4[1][1][2][2] + 
              2*xz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[1][1][2][2] + 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[1][1][2][2] - 
              2*xx*pow(param->gamm[1*3+2],2)*a4[1][1][2][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[1][1][2][2] - 
              2*xz*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[1][1][2][2] - 
              2*xx*param->gamm[1*3+1]*param->gamm[2*3+2]*a4[1][1][2][2] + 
              xx*pow(param->gamm[2*3+2],2)*a4[1][1][2][2] + 
              2*yy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[1][2][2][2] - 
              2*yz*pow(param->gamm[0*3+2],2)*a4[1][2][2][2] - 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[1][2][2][2] - 
              2*xy*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[1][2][2][2] + 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[1][2][2][2] + 
              2*xx*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[1][2][2][2] + 
              2*xy*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[1][2][2][2] - 
              2*xx*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[1][2][2][2] + 
              yy*pow(param->gamm[0*3+2],2)*a4[2][2][2][2] - 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[2][2][2][2] + 
              xx*pow(param->gamm[1*3+2],2)*a4[2][2][2][2];
      psidot[ind(l,m,0)] += -l*(l+1)*(param->C1/param->normgamma*temp1/4 + param->C2*param->normgamma*psi[ind(l,m,0)]);
      xx = (mc[0]*psi[ind(l-2,m-2,1)]+mc[1]*psi[ind(l-2,m,1)]+mc[2]*psi[ind(l-2,m+2,1)]+mc[3]*psi[ind(l,m-2,1)]+mc[4]*psi[ind(l,m,1)]+mc[5]*psi[ind(l,m+2,1)]+mc[6]*psi[ind(l+2,m-2,1)]+mc[7]*psi[ind(l+2,m,1)]+mc[8]*psi[ind(l+2,m+2,1)]);
      xy = (mc[9]*psi[ind(l-2,m-2,0)]+mc[10]*psi[ind(l-2,m+2,0)]+mc[11]*psi[ind(l,m-2,0)]+mc[12]*psi[ind(l,m+2,0)]+mc[13]*psi[ind(l+2,m-2,0)]+mc[14]*psi[ind(l+2,m+2,0)]);
      xz = (mc[15]*psi[ind(l-2,m-1,1)]+mc[16]*psi[ind(l-2,m+1,1)]+mc[17]*psi[ind(l,m-1,1)]+mc[18]*psi[ind(l,m+1,1)]+mc[19]*psi[ind(l+2,m-1,1)]+mc[20]*psi[ind(l+2,m+1,1)]);
      yy = (mc[21]*psi[ind(l-2,m-2,1)]+mc[22]*psi[ind(l-2,m,1)]+mc[23]*psi[ind(l-2,m+2,1)]+mc[24]*psi[ind(l,m-2,1)]+mc[25]*psi[ind(l,m,1)]+mc[26]*psi[ind(l,m+2,1)]+mc[27]*psi[ind(l+2,m-2,1)]+mc[28]*psi[ind(l+2,m,1)]+mc[29]*psi[ind(l+2,m+2,1)]);
      yz = (mc[30]*psi[ind(l-2,m-1,0)]+mc[31]*psi[ind(l-2,m+1,0)]+mc[32]*psi[ind(l,m-1,0)]+mc[33]*psi[ind(l,m+1,0)]+mc[34]*psi[ind(l+2,m-1,0)]+mc[35]*psi[ind(l+2,m+1,0)]);
      zz = (mc[36]*psi[ind(l-2,m,1)]+mc[37]*psi[ind(l,m,1)]+mc[38]*psi[ind(l+2,m,1)]);
      xxxx = (mc[39]*psi[ind(l-2,m-2,1)]+mc[40]*psi[ind(l-2,m-4,1)]+mc[41]*psi[ind(l-2,m,1)]+mc[42]*psi[ind(l-2,m+2,1)]+mc[43]*psi[ind(l-2,m+4,1)]+mc[44]*psi[ind(l-4,m-2,1)]+mc[45]*psi[ind(l-4,m-4,1)]+mc[46]*psi[ind(l-4,m,1)]+mc[47]*psi[ind(l-4,m+2,1)]+mc[48]*psi[ind(l-4,m+4,1)]+mc[49]*psi[ind(l,m-2,1)]+mc[50]*psi[ind(l,m-4,1)]+mc[51]*psi[ind(l,m,1)]+mc[52]*psi[ind(l,m+2,1)]+mc[53]*psi[ind(l,m+4,1)]+mc[54]*psi[ind(l+2,m-2,1)]+mc[55]*psi[ind(l+2,m-4,1)]+mc[56]*psi[ind(l+2,m,1)]+mc[57]*psi[ind(l+2,m+2,1)]+mc[58]*psi[ind(l+2,m+4,1)]+mc[59]*psi[ind(l+4,m-2,1)]+mc[60]*psi[ind(l+4,m-4,1)]+mc[61]*psi[ind(l+4,m,1)]+mc[62]*psi[ind(l+4,m+2,1)]+mc[63]*psi[ind(l+4,m+4,1)]);
      xxxy = (mc[64]*psi[ind(l-2,m-2,0)]+mc[65]*psi[ind(l-2,m-4,0)]+mc[66]*psi[ind(l-2,m+2,0)]+mc[67]*psi[ind(l-2,m+4,0)]+mc[68]*psi[ind(l-4,m-2,0)]+mc[69]*psi[ind(l-4,m-4,0)]+mc[70]*psi[ind(l-4,m+2,0)]+mc[71]*psi[ind(l-4,m+4,0)]+mc[72]*psi[ind(l,m-2,0)]+mc[73]*psi[ind(l,m-4,0)]+mc[74]*psi[ind(l,m+2,0)]+mc[75]*psi[ind(l,m+4,0)]+mc[76]*psi[ind(l+2,m-2,0)]+mc[77]*psi[ind(l+2,m-4,0)]+mc[78]*psi[ind(l+2,m+2,0)]+mc[79]*psi[ind(l+2,m+4,0)]+mc[80]*psi[ind(l+4,m-2,0)]+mc[81]*psi[ind(l+4,m-4,0)]+mc[82]*psi[ind(l+4,m+2,0)]+mc[83]*psi[ind(l+4,m+4,0)]);
      xxxz = (mc[84]*psi[ind(l-2,m-1,1)]+mc[85]*psi[ind(l-2,m-3,1)]+mc[86]*psi[ind(l-2,m+1,1)]+mc[87]*psi[ind(l-2,m+3,1)]+mc[88]*psi[ind(l-4,m-1,1)]+mc[89]*psi[ind(l-4,m-3,1)]+mc[90]*psi[ind(l-4,m+1,1)]+mc[91]*psi[ind(l-4,m+3,1)]+mc[92]*psi[ind(l,m-1,1)]+mc[93]*psi[ind(l,m-3,1)]+mc[94]*psi[ind(l,m+1,1)]+mc[95]*psi[ind(l,m+3,1)]+mc[96]*psi[ind(l+2,m-1,1)]+mc[97]*psi[ind(l+2,m-3,1)]+mc[98]*psi[ind(l+2,m+1,1)]+mc[99]*psi[ind(l+2,m+3,1)]+mc[100]*psi[ind(l+4,m-1,1)]+mc[101]*psi[ind(l+4,m-3,1)]+mc[102]*psi[ind(l+4,m+1,1)]+mc[103]*psi[ind(l+4,m+3,1)]);
      xxyy = (mc[104]*psi[ind(l-2,m-4,1)]+mc[105]*psi[ind(l-2,m,1)]+mc[106]*psi[ind(l-2,m+4,1)]+mc[107]*psi[ind(l-4,m-4,1)]+mc[108]*psi[ind(l-4,m,1)]+mc[109]*psi[ind(l-4,m+4,1)]+mc[110]*psi[ind(l,m-4,1)]+mc[111]*psi[ind(l,m,1)]+mc[112]*psi[ind(l,m+4,1)]+mc[113]*psi[ind(l+2,m-4,1)]+mc[114]*psi[ind(l+2,m,1)]+mc[115]*psi[ind(l+2,m+4,1)]+mc[116]*psi[ind(l+4,m-4,1)]+mc[117]*psi[ind(l+4,m,1)]+mc[118]*psi[ind(l+4,m+4,1)]);
      xxyz = (mc[119]*psi[ind(l-2,m-1,0)]+mc[120]*psi[ind(l-2,m-3,0)]+mc[121]*psi[ind(l-2,m+1,0)]+mc[122]*psi[ind(l-2,m+3,0)]+mc[123]*psi[ind(l-4,m-1,0)]+mc[124]*psi[ind(l-4,m-3,0)]+mc[125]*psi[ind(l-4,m+1,0)]+mc[126]*psi[ind(l-4,m+3,0)]+mc[127]*psi[ind(l,m-1,0)]+mc[128]*psi[ind(l,m-3,0)]+mc[129]*psi[ind(l,m+1,0)]+mc[130]*psi[ind(l,m+3,0)]+mc[131]*psi[ind(l+2,m-1,0)]+mc[132]*psi[ind(l+2,m-3,0)]+mc[133]*psi[ind(l+2,m+1,0)]+mc[134]*psi[ind(l+2,m+3,0)]+mc[135]*psi[ind(l+4,m-1,0)]+mc[136]*psi[ind(l+4,m-3,0)]+mc[137]*psi[ind(l+4,m+1,0)]+mc[138]*psi[ind(l+4,m+3,0)]);
      xxzz = (mc[139]*psi[ind(l-2,m-2,1)]+mc[140]*psi[ind(l-2,m,1)]+mc[141]*psi[ind(l-2,m+2,1)]+mc[142]*psi[ind(l-4,m-2,1)]+mc[143]*psi[ind(l-4,m,1)]+mc[144]*psi[ind(l-4,m+2,1)]+mc[145]*psi[ind(l,m-2,1)]+mc[146]*psi[ind(l,m,1)]+mc[147]*psi[ind(l,m+2,1)]+mc[148]*psi[ind(l+2,m-2,1)]+mc[149]*psi[ind(l+2,m,1)]+mc[150]*psi[ind(l+2,m+2,1)]+mc[151]*psi[ind(l+4,m-2,1)]+mc[152]*psi[ind(l+4,m,1)]+mc[153]*psi[ind(l+4,m+2,1)]);
      xyyy = (mc[154]*psi[ind(l-2,m-2,0)]+mc[155]*psi[ind(l-2,m-4,0)]+mc[156]*psi[ind(l-2,m+2,0)]+mc[157]*psi[ind(l-2,m+4,0)]+mc[158]*psi[ind(l-4,m-2,0)]+mc[159]*psi[ind(l-4,m-4,0)]+mc[160]*psi[ind(l-4,m+2,0)]+mc[161]*psi[ind(l-4,m+4,0)]+mc[162]*psi[ind(l,m-2,0)]+mc[163]*psi[ind(l,m-4,0)]+mc[164]*psi[ind(l,m+2,0)]+mc[165]*psi[ind(l,m+4,0)]+mc[166]*psi[ind(l+2,m-2,0)]+mc[167]*psi[ind(l+2,m-4,0)]+mc[168]*psi[ind(l+2,m+2,0)]+mc[169]*psi[ind(l+2,m+4,0)]+mc[170]*psi[ind(l+4,m-2,0)]+mc[171]*psi[ind(l+4,m-4,0)]+mc[172]*psi[ind(l+4,m+2,0)]+mc[173]*psi[ind(l+4,m+4,0)]);
      xyyz = (mc[174]*psi[ind(l-2,m-1,1)]+mc[175]*psi[ind(l-2,m-3,1)]+mc[176]*psi[ind(l-2,m+1,1)]+mc[177]*psi[ind(l-2,m+3,1)]+mc[178]*psi[ind(l-4,m-1,1)]+mc[179]*psi[ind(l-4,m-3,1)]+mc[180]*psi[ind(l-4,m+1,1)]+mc[181]*psi[ind(l-4,m+3,1)]+mc[182]*psi[ind(l,m-1,1)]+mc[183]*psi[ind(l,m-3,1)]+mc[184]*psi[ind(l,m+1,1)]+mc[185]*psi[ind(l,m+3,1)]+mc[186]*psi[ind(l+2,m-1,1)]+mc[187]*psi[ind(l+2,m-3,1)]+mc[188]*psi[ind(l+2,m+1,1)]+mc[189]*psi[ind(l+2,m+3,1)]+mc[190]*psi[ind(l+4,m-1,1)]+mc[191]*psi[ind(l+4,m-3,1)]+mc[192]*psi[ind(l+4,m+1,1)]+mc[193]*psi[ind(l+4,m+3,1)]);
      xyzz = (mc[194]*psi[ind(l-2,m-2,0)]+mc[195]*psi[ind(l-2,m+2,0)]+mc[196]*psi[ind(l-4,m-2,0)]+mc[197]*psi[ind(l-4,m+2,0)]+mc[198]*psi[ind(l,m-2,0)]+mc[199]*psi[ind(l,m+2,0)]+mc[200]*psi[ind(l+2,m-2,0)]+mc[201]*psi[ind(l+2,m+2,0)]+mc[202]*psi[ind(l+4,m-2,0)]+mc[203]*psi[ind(l+4,m+2,0)]);
      xzzz = (mc[204]*psi[ind(l-2,m-1,1)]+mc[205]*psi[ind(l-2,m+1,1)]+mc[206]*psi[ind(l-4,m-1,1)]+mc[207]*psi[ind(l-4,m+1,1)]+mc[208]*psi[ind(l,m-1,1)]+mc[209]*psi[ind(l,m+1,1)]+mc[210]*psi[ind(l+2,m-1,1)]+mc[211]*psi[ind(l+2,m+1,1)]+mc[212]*psi[ind(l+4,m-1,1)]+mc[213]*psi[ind(l+4,m+1,1)]);
      yyyy = (mc[214]*psi[ind(l-2,m-2,1)]+mc[215]*psi[ind(l-2,m-4,1)]+mc[216]*psi[ind(l-2,m,1)]+mc[217]*psi[ind(l-2,m+2,1)]+mc[218]*psi[ind(l-2,m+4,1)]+mc[219]*psi[ind(l-4,m-2,1)]+mc[220]*psi[ind(l-4,m-4,1)]+mc[221]*psi[ind(l-4,m,1)]+mc[222]*psi[ind(l-4,m+2,1)]+mc[223]*psi[ind(l-4,m+4,1)]+mc[224]*psi[ind(l,m-2,1)]+mc[225]*psi[ind(l,m-4,1)]+mc[226]*psi[ind(l,m,1)]+mc[227]*psi[ind(l,m+2,1)]+mc[228]*psi[ind(l,m+4,1)]+mc[229]*psi[ind(l+2,m-2,1)]+mc[230]*psi[ind(l+2,m-4,1)]+mc[231]*psi[ind(l+2,m,1)]+mc[232]*psi[ind(l+2,m+2,1)]+mc[233]*psi[ind(l+2,m+4,1)]+mc[234]*psi[ind(l+4,m-2,1)]+mc[235]*psi[ind(l+4,m-4,1)]+mc[236]*psi[ind(l+4,m,1)]+mc[237]*psi[ind(l+4,m+2,1)]+mc[238]*psi[ind(l+4,m+4,1)]);
      yyyz = (mc[239]*psi[ind(l-2,m-1,0)]+mc[240]*psi[ind(l-2,m-3,0)]+mc[241]*psi[ind(l-2,m+1,0)]+mc[242]*psi[ind(l-2,m+3,0)]+mc[243]*psi[ind(l-4,m-1,0)]+mc[244]*psi[ind(l-4,m-3,0)]+mc[245]*psi[ind(l-4,m+1,0)]+mc[246]*psi[ind(l-4,m+3,0)]+mc[247]*psi[ind(l,m-1,0)]+mc[248]*psi[ind(l,m-3,0)]+mc[249]*psi[ind(l,m+1,0)]+mc[250]*psi[ind(l,m+3,0)]+mc[251]*psi[ind(l+2,m-1,0)]+mc[252]*psi[ind(l+2,m-3,0)]+mc[253]*psi[ind(l+2,m+1,0)]+mc[254]*psi[ind(l+2,m+3,0)]+mc[255]*psi[ind(l+4,m-1,0)]+mc[256]*psi[ind(l+4,m-3,0)]+mc[257]*psi[ind(l+4,m+1,0)]+mc[258]*psi[ind(l+4,m+3,0)]);
      yyzz = (mc[259]*psi[ind(l-2,m-2,1)]+mc[260]*psi[ind(l-2,m,1)]+mc[261]*psi[ind(l-2,m+2,1)]+mc[262]*psi[ind(l-4,m-2,1)]+mc[263]*psi[ind(l-4,m,1)]+mc[264]*psi[ind(l-4,m+2,1)]+mc[265]*psi[ind(l,m-2,1)]+mc[266]*psi[ind(l,m,1)]+mc[267]*psi[ind(l,m+2,1)]+mc[268]*psi[ind(l+2,m-2,1)]+mc[269]*psi[ind(l+2,m,1)]+mc[270]*psi[ind(l+2,m+2,1)]+mc[271]*psi[ind(l+4,m-2,1)]+mc[272]*psi[ind(l+4,m,1)]+mc[273]*psi[ind(l+4,m+2,1)]);
      yzzz = (mc[274]*psi[ind(l-2,m-1,0)]+mc[275]*psi[ind(l-2,m+1,0)]+mc[276]*psi[ind(l-4,m-1,0)]+mc[277]*psi[ind(l-4,m+1,0)]+mc[278]*psi[ind(l,m-1,0)]+mc[279]*psi[ind(l,m+1,0)]+mc[280]*psi[ind(l+2,m-1,0)]+mc[281]*psi[ind(l+2,m+1,0)]+mc[282]*psi[ind(l+4,m-1,0)]+mc[283]*psi[ind(l+4,m+1,0)]);
      zzzz = (mc[284]*psi[ind(l-2,m,1)]+mc[285]*psi[ind(l-4,m,1)]+mc[286]*psi[ind(l,m,1)]+mc[287]*psi[ind(l+2,m,1)]+mc[288]*psi[ind(l+4,m,1)]);
      temp1 = xxzz*pow(param->gamm[0*3+1],2)*a2[0][0] - 
              2*xxyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][0] + 
              xxyy*pow(param->gamm[0*3+2],2)*a2[0][0] + 
              2*xyzz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[0][0] - 
              2*xyyz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][0] + 
              yyzz*pow(param->gamm[1*3+1],2)*a2[0][0] - 
              2*xyyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][0] + 
              2*xzzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][0] + 
              2*xyyy*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][0] - 
              2*xyzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][0] - 
              2*yyyz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][0] + 
              2*yzzz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][0] + 
              yyyy*pow(param->gamm[1*3+2],2)*a2[0][0] - 
              2*yyzz*pow(param->gamm[1*3+2],2)*a2[0][0] + 
              zzzz*pow(param->gamm[1*3+2],2)*a2[0][0] - 
              2*xyzz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][0] + 
              2*xyyz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[0][0] - 
              2*yyzz*param->gamm[1*3+1]*param->gamm[2*3+2]*a2[0][0] + 
              2*yyyz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][0] - 
              2*yzzz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][0] + 
              yyzz*pow(param->gamm[2*3+2],2)*a2[0][0] - 
              2*xxzz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[0][1] - 
              2*xyzz*pow(param->gamm[0*3+1],2)*a2[0][1] + 
              2*xxyz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[0][1] + 
              2*xxxz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][1] + 
              2*xyyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][1] - 
              2*xzzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][1] - 
              2*xxxy*pow(param->gamm[0*3+2],2)*a2[0][1] + 
              2*xyzz*pow(param->gamm[0*3+2],2)*a2[0][1] - 
              2*xyzz*param->gamm[0*3+0]*param->gamm[1*3+1]*a2[0][1] - 
              2*yyzz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[0][1] + 
              2*xxyz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][1] - 
              2*yzzz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][1] + 
              2*xyyz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[0][1] - 
              2*xzzz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[0][1] + 
              2*xxyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][1] + 
              2*yyyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][1] - 
              2*yzzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][1] - 
              4*xxyy*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][1] + 
              2*xxzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][1] + 
              2*yyzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][1] - 
              2*zzzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][1] + 
              2*xyyz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][1] - 
              2*xyyy*pow(param->gamm[1*3+2],2)*a2[0][1] + 
              2*xyzz*pow(param->gamm[1*3+2],2)*a2[0][1] + 
              2*xyzz*param->gamm[0*3+0]*param->gamm[2*3+2]*a2[0][1] + 
              2*xxzz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][1] + 
              2*yyzz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][1] - 
              4*xxyz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[0][1] + 
              2*yzzz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[0][1] + 
              2*xyzz*param->gamm[1*3+1]*param->gamm[2*3+2]*a2[0][1] - 
              4*xyyz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][1] + 
              2*xzzz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][1] - 
              2*xyzz*pow(param->gamm[2*3+2],2)*a2[0][1] + 
              2*xxyz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[0][2] - 
              2*xxxz*pow(param->gamm[0*3+1],2)*a2[0][2] + 
              2*xyyz*pow(param->gamm[0*3+1],2)*a2[0][2] - 
              2*xxyy*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[0][2] + 
              2*xxxy*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][2] - 
              2*xyyy*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][2] + 
              2*xyzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[0][2] - 
              2*xyyz*pow(param->gamm[0*3+2],2)*a2[0][2] + 
              2*xyyz*param->gamm[0*3+0]*param->gamm[1*3+1]*a2[0][2] - 
              4*xxyz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[0][2] + 
              2*yyyz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[0][2] + 
              2*xxyy*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][2] + 
              2*yyzz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[0][2] - 
              2*xyyz*pow(param->gamm[1*3+1],2)*a2[0][2] - 
              2*xyyy*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[0][2] + 
              2*xyzz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[0][2] + 
              2*xxyy*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][2] - 
              4*xxzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][2] - 
              2*yyyy*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][2] + 
              2*yyzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[0][2] + 
              2*xxyz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][2] - 
              2*yyyz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][2] + 
              2*yzzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[0][2] + 
              2*xyyy*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][2] - 
              4*xyzz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[0][2] + 
              2*xyyz*pow(param->gamm[1*3+2],2)*a2[0][2] - 
              2*xzzz*pow(param->gamm[1*3+2],2)*a2[0][2] - 
              2*xyyz*param->gamm[0*3+0]*param->gamm[2*3+2]*a2[0][2] + 
              2*xxyz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][2] - 
              2*yyyz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[0][2] - 
              2*yyzz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[0][2] + 
              2*xyyz*param->gamm[1*3+1]*param->gamm[2*3+2]*a2[0][2] + 
              2*xyzz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[0][2] + 
              xxzz*pow(param->gamm[0*3+0],2)*a2[1][1] + 
              2*xyzz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[1][1] + 
              yyzz*pow(param->gamm[0*3+1],2)*a2[1][1] - 
              2*xxxz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[1][1] + 
              2*xzzz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[1][1] - 
              2*xxyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][1] + 
              2*yzzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][1] + 
              xxxx*pow(param->gamm[0*3+2],2)*a2[1][1] - 
              2*xxzz*pow(param->gamm[0*3+2],2)*a2[1][1] + 
              zzzz*pow(param->gamm[0*3+2],2)*a2[1][1] - 
              2*xxyz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[1][1] - 
              2*xyyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[1][1] + 
              2*xxxy*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][1] - 
              2*xyzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][1] + 
              xxyy*pow(param->gamm[1*3+2],2)*a2[1][1] - 
              2*xxzz*param->gamm[0*3+0]*param->gamm[2*3+2]*a2[1][1] - 
              2*xyzz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[1][1] + 
              2*xxxz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[1][1] - 
              2*xzzz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[1][1] + 
              2*xxyz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[1][1] + 
              xxzz*pow(param->gamm[2*3+2],2)*a2[1][1] - 
              2*xxyz*pow(param->gamm[0*3+0],2)*a2[1][2] + 
              2*xxxz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[1][2] - 
              4*xyyz*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[1][2] + 
              2*xxyz*pow(param->gamm[0*3+1],2)*a2[1][2] - 
              2*yyyz*pow(param->gamm[0*3+1],2)*a2[1][2] + 
              2*xxxy*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[1][2] - 
              4*xyzz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[1][2] - 
              2*xxxx*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][2] + 
              2*xxyy*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][2] + 
              2*xxzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][2] - 
              4*yyzz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[1][2] + 
              2*xxyz*pow(param->gamm[0*3+2],2)*a2[1][2] - 
              2*yzzz*pow(param->gamm[0*3+2],2)*a2[1][2] + 
              2*xxyz*param->gamm[0*3+0]*param->gamm[1*3+1]*a2[1][2] + 
              2*xyyz*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[1][2] - 
              2*xxxy*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[1][2] + 
              2*xyzz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[1][2] + 
              2*xxyy*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[1][2] + 
              2*xxzz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[1][2] - 
              2*xxxy*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[1][2] + 
              2*xyyy*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[1][2] + 
              2*xyzz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[1][2] - 
              2*xxxz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][2] + 
              2*xyyz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][2] + 
              2*xzzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[1][2] - 
              2*xxyy*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[1][2] - 
              2*xxyz*pow(param->gamm[1*3+2],2)*a2[1][2] + 
              2*xxyz*param->gamm[0*3+0]*param->gamm[2*3+2]*a2[1][2] - 
              2*xxxz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[1][2] + 
              2*xyyz*param->gamm[0*3+1]*param->gamm[2*3+2]*a2[1][2] + 
              2*xyzz*param->gamm[0*3+2]*param->gamm[2*3+2]*a2[1][2] - 
              2*xxyz*param->gamm[1*3+1]*param->gamm[2*3+2]*a2[1][2] - 
              2*xxzz*param->gamm[1*3+2]*param->gamm[2*3+2]*a2[1][2] + 
              xxyy*pow(param->gamm[0*3+0],2)*a2[2][2] - 
              2*xxxy*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[2][2] + 
              2*xyyy*param->gamm[0*3+0]*param->gamm[0*3+1]*a2[2][2] + 
              xxxx*pow(param->gamm[0*3+1],2)*a2[2][2] - 
              2*xxyy*pow(param->gamm[0*3+1],2)*a2[2][2] + 
              yyyy*pow(param->gamm[0*3+1],2)*a2[2][2] + 
              2*xyyz*param->gamm[0*3+0]*param->gamm[0*3+2]*a2[2][2] - 
              2*xxyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[2][2] + 
              2*yyyz*param->gamm[0*3+1]*param->gamm[0*3+2]*a2[2][2] + 
              yyzz*pow(param->gamm[0*3+2],2)*a2[2][2] - 
              2*xxyy*param->gamm[0*3+0]*param->gamm[1*3+1]*a2[2][2] + 
              2*xxxy*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[2][2] - 
              2*xyyy*param->gamm[0*3+1]*param->gamm[1*3+1]*a2[2][2] - 
              2*xyyz*param->gamm[0*3+2]*param->gamm[1*3+1]*a2[2][2] + 
              xxyy*pow(param->gamm[1*3+1],2)*a2[2][2] - 
              2*xxyz*param->gamm[0*3+0]*param->gamm[1*3+2]*a2[2][2] + 
              2*xxxz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[2][2] - 
              2*xyyz*param->gamm[0*3+1]*param->gamm[1*3+2]*a2[2][2] - 
              2*xyzz*param->gamm[0*3+2]*param->gamm[1*3+2]*a2[2][2] + 
              2*xxyz*param->gamm[1*3+1]*param->gamm[1*3+2]*a2[2][2] + 
              xxzz*pow(param->gamm[1*3+2],2)*a2[2][2] + 
              zz*pow(param->gamm[0*3+1],2)*a4[0][0][0][0] - 
              2*yz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][0][0] + 
              yy*pow(param->gamm[0*3+2],2)*a4[0][0][0][0] - 
              2*zz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][0][0][1] + 
              2*yz*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][0][0][1] + 
              2*xz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][0][1] - 
              2*xy*pow(param->gamm[0*3+2],2)*a4[0][0][0][1] + 
              2*zz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][0][0][1] - 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][0][0][1] - 
              2*yz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][0][1] + 
              2*yy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][0][1] + 
              2*yz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][0][0][2] - 
              2*xz*pow(param->gamm[0*3+1],2)*a4[0][0][0][2] - 
              2*yy*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][0][0][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][0][2] + 
              2*zz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][0][2] - 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][0][2] - 
              2*yz*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][0][0][2] + 
              2*yy*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][0][0][2] + 
              zz*pow(param->gamm[0*3+0],2)*a4[0][0][1][1] - 
              2*zz*pow(param->gamm[0*3+1],2)*a4[0][0][1][1] - 
              2*xz*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][0][1][1] + 
              2*yz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][1][1] + 
              xx*pow(param->gamm[0*3+2],2)*a4[0][0][1][1] - 
              2*zz*param->gamm[0*3+0]*param->gamm[1*3+1]*a4[0][0][1][1] + 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][0][1][1] + 
              zz*pow(param->gamm[1*3+1],2)*a4[0][0][1][1] + 
              2*yz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][0][1][1] + 
              2*xz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][1][1] - 
              4*xy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][1][1] - 
              2*yz*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][0][1][1] + 
              yy*pow(param->gamm[1*3+2],2)*a4[0][0][1][1] - 
              2*yz*pow(param->gamm[0*3+0],2)*a4[0][0][1][2] + 
              2*xz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][0][1][2] + 
              2*yz*pow(param->gamm[0*3+1],2)*a4[0][0][1][2] + 
              2*xy*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][0][1][2] - 
              2*xx*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][1][2] - 
              2*yy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][1][2] - 
              2*zz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][1][2] + 
              2*yz*pow(param->gamm[0*3+2],2)*a4[0][0][1][2] + 
              2*yz*param->gamm[0*3+0]*param->gamm[1*3+1]*a4[0][0][1][2] - 
              4*xz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][0][1][2] + 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][0][1][2] - 
              2*yy*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][0][1][2] - 
              2*zz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][0][1][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][1][2] + 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][1][2] + 
              2*zz*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][0][1][2] - 
              2*yz*pow(param->gamm[1*3+2],2)*a4[0][0][1][2] + 
              2*yz*param->gamm[0*3+0]*param->gamm[2*3+2]*a4[0][0][1][2] + 
              2*xz*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][0][1][2] - 
              4*xy*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][0][1][2] - 
              2*yz*param->gamm[1*3+1]*param->gamm[2*3+2]*a4[0][0][1][2] + 
              2*yy*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][0][1][2] + 
              yy*pow(param->gamm[0*3+0],2)*a4[0][0][2][2] - 
              2*xy*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][0][2][2] + 
              xx*pow(param->gamm[0*3+1],2)*a4[0][0][2][2] + 
              2*yz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][0][2][2] - 
              2*yy*pow(param->gamm[0*3+2],2)*a4[0][0][2][2] + 
              2*yz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][0][2][2] - 
              4*xz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][0][2][2] + 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][0][2][2] + 
              zz*pow(param->gamm[1*3+2],2)*a4[0][0][2][2] - 
              2*yy*param->gamm[0*3+0]*param->gamm[2*3+2]*a4[0][0][2][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][0][2][2] - 
              2*yz*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][0][2][2] + 
              yy*pow(param->gamm[2*3+2],2)*a4[0][0][2][2] + 
              2*zz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][1][1][1] - 
              2*xz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][1][1][1] - 
              2*zz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][1][1][1] - 
              2*xz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][1][1][1] + 
              2*yz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][1][1] + 
              2*xx*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][1][1] + 
              2*xz*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][1][1][1] - 
              2*xy*pow(param->gamm[1*3+2],2)*a4[0][1][1][1] - 
              4*yz*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][1][1][2] + 
              2*xz*pow(param->gamm[0*3+1],2)*a4[0][1][1][2] + 
              2*zz*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][1][1][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][1][1][2] - 
              2*xz*pow(param->gamm[0*3+2],2)*a4[0][1][1][2] + 
              2*xz*param->gamm[0*3+0]*param->gamm[1*3+1]*a4[0][1][1][2] + 
              2*yz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][1][1][2] - 
              2*xx*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][1][1][2] - 
              2*zz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][1][1][2] - 
              2*xz*pow(param->gamm[1*3+1],2)*a4[0][1][1][2] + 
              2*xy*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][1][1][2] - 
              2*xx*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][1][2] - 
              2*yy*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][1][2] - 
              2*zz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][1][2] + 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][1][2] + 
              2*xy*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][1][1][2] + 
              2*xz*pow(param->gamm[1*3+2],2)*a4[0][1][1][2] - 
              2*xz*param->gamm[0*3+0]*param->gamm[2*3+2]*a4[0][1][1][2] + 
              2*yz*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][1][1][2] + 
              2*xx*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][1][1][2] + 
              2*xz*param->gamm[1*3+1]*param->gamm[2*3+2]*a4[0][1][1][2] - 
              4*xy*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][1][1][2] + 
              2*yy*param->gamm[0*3+0]*param->gamm[0*3+1]*a4[0][1][2][2] - 
              2*xy*pow(param->gamm[0*3+1],2)*a4[0][1][2][2] - 
              4*yz*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][1][2][2] + 
              2*xz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][1][2][2] + 
              2*xy*pow(param->gamm[0*3+2],2)*a4[0][1][2][2] - 
              2*xy*param->gamm[0*3+0]*param->gamm[1*3+1]*a4[0][1][2][2] + 
              2*xx*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[0][1][2][2] + 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[0][1][2][2] + 
              2*xz*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][1][2][2] + 
              2*yz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][1][2][2] - 
              2*xx*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][2][2] - 
              2*yy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][2][2] - 
              2*zz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][1][2][2] - 
              4*xz*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[0][1][2][2] + 
              2*xy*pow(param->gamm[1*3+2],2)*a4[0][1][2][2] + 
              2*xy*param->gamm[0*3+0]*param->gamm[2*3+2]*a4[0][1][2][2] - 
              2*xx*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][1][2][2] - 
              2*yy*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[0][1][2][2] + 
              2*yz*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][1][2][2] + 
              2*xy*param->gamm[1*3+1]*param->gamm[2*3+2]*a4[0][1][2][2] + 
              2*xz*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][1][2][2] - 
              2*xy*pow(param->gamm[2*3+2],2)*a4[0][1][2][2] + 
              2*yy*param->gamm[0*3+0]*param->gamm[0*3+2]*a4[0][2][2][2] - 
              2*xy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[0][2][2][2] - 
              2*xy*param->gamm[0*3+0]*param->gamm[1*3+2]*a4[0][2][2][2] + 
              2*xx*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[0][2][2][2] + 
              2*yz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[0][2][2][2] - 
              2*xz*pow(param->gamm[1*3+2],2)*a4[0][2][2][2] - 
              2*yy*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[0][2][2][2] + 
              2*xy*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[0][2][2][2] + 
              zz*pow(param->gamm[0*3+1],2)*a4[1][1][1][1] - 
              2*xz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[1][1][1][1] + 
              xx*pow(param->gamm[1*3+2],2)*a4[1][1][1][1] - 
              2*yz*pow(param->gamm[0*3+1],2)*a4[1][1][1][2] + 
              2*zz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[1][1][1][2] + 
              2*xz*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[1][1][1][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[1][1][1][2] - 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[1][1][1][2] - 
              2*xx*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[1][1][1][2] - 
              2*xz*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[1][1][1][2] + 
              2*xx*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[1][1][1][2] + 
              yy*pow(param->gamm[0*3+1],2)*a4[1][1][2][2] - 
              4*yz*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[1][1][2][2] + 
              zz*pow(param->gamm[0*3+2],2)*a4[1][1][2][2] - 
              2*xy*param->gamm[0*3+1]*param->gamm[1*3+1]*a4[1][1][2][2] + 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[1][1][2][2] + 
              xx*pow(param->gamm[1*3+1],2)*a4[1][1][2][2] + 
              2*xz*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[1][1][2][2] + 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[1][1][2][2] - 
              2*xx*pow(param->gamm[1*3+2],2)*a4[1][1][2][2] + 
              2*xy*param->gamm[0*3+1]*param->gamm[2*3+2]*a4[1][1][2][2] - 
              2*xz*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[1][1][2][2] - 
              2*xx*param->gamm[1*3+1]*param->gamm[2*3+2]*a4[1][1][2][2] + 
              xx*pow(param->gamm[2*3+2],2)*a4[1][1][2][2] + 
              2*yy*param->gamm[0*3+1]*param->gamm[0*3+2]*a4[1][2][2][2] - 
              2*yz*pow(param->gamm[0*3+2],2)*a4[1][2][2][2] - 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+1]*a4[1][2][2][2] - 
              2*xy*param->gamm[0*3+1]*param->gamm[1*3+2]*a4[1][2][2][2] + 
              2*xz*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[1][2][2][2] + 
              2*xx*param->gamm[1*3+1]*param->gamm[1*3+2]*a4[1][2][2][2] + 
              2*xy*param->gamm[0*3+2]*param->gamm[2*3+2]*a4[1][2][2][2] - 
              2*xx*param->gamm[1*3+2]*param->gamm[2*3+2]*a4[1][2][2][2] + 
              yy*pow(param->gamm[0*3+2],2)*a4[2][2][2][2] - 
              2*xy*param->gamm[0*3+2]*param->gamm[1*3+2]*a4[2][2][2][2] + 
              xx*pow(param->gamm[1*3+2],2)*a4[2][2][2][2];
      psidot[ind(l,m,1)] += -l*(l+1)*(param->C1/param->normgamma*temp1/4 + param->C2*param->normgamma*psi[ind(l,m,1)]);
    }
    /* Broadcast that thread is finished. */
    sem_post(job_info->sem_end);
  }
  return(0);
}
