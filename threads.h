#ifdef WIN32

#include <windows.h>
#define sem_t HANDLE
#define sem_init(sem) sem = CreateSemaphore(NULL, 0, 256000, NULL)
#define sem_post(sem) ReleaseSemaphore(sem,1,NULL)
#define sem_wait(sem) WaitForSingleObject(sem, INFINITE)
static DWORD tid;
#define thread_create(function,data) CreateThread(NULL,0,function,data,0,&tid);
#define thread_return_t DWORD WINAPI

#else

#include <pthread.h>
static pthread_mutex_t job_m = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t job_c = PTHREAD_COND_INITIALIZER;
#define sem_t int
#define sem_init(sem) sem = 0;
#define sem_post(sem)                           \
  do {                                          \
    pthread_mutex_lock(&job_m);                 \
    sem+=1;                                     \
    pthread_mutex_unlock(&job_m);               \
    pthread_cond_broadcast(&job_c);             \
  } while(0)
#define sem_wait(sem)                           \
  do {                                          \
    pthread_mutex_lock(&job_m);                 \
    while(sem==0)                               \
      pthread_cond_wait(&job_c, &job_m);        \
    sem--;                                      \
    pthread_mutex_unlock(&job_m);               \
  } while(0)
static pthread_t pid;
#define thread_create(function,data)            \
  do {                                          \
    pthread_create(&pid, NULL, function, data); \
    pthread_detach(pid);                        \
  } while(0)
#define thread_return_t void *

#endif

