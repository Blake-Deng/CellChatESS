#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>

/* Observe real GNU OpenMP teams without changing ESS source or its DLL.
 * Used only in correctness validation, never in performance measurements. */
typedef void (*parallel_fn)(void (*)(void *), void *, unsigned, unsigned);
struct probe_data { void (*fn)(void *); void *data; unsigned requested; const char *caller; };

static void probe_worker(void *ptr) {
  struct probe_data *p=(struct probe_data *)ptr;
  int (*team_size)(void)=dlsym(RTLD_DEFAULT,"omp_get_num_threads");
  int (*thread_id)(void)=dlsym(RTLD_DEFAULT,"omp_get_thread_num");
  const char *path=getenv("CELLCHATESS_OPENMP_TEAM_LOG");
  if(path && team_size && thread_id && thread_id()==0) {
    char line[4096];
    int n=snprintf(line,sizeof(line),"%d\t%u\t%d\t%s\n",getpid(),p->requested,team_size(),p->caller);
    int fd=open(path,O_WRONLY|O_CREAT|O_APPEND,0600);
    if(fd>=0) {write(fd,line,n);close(fd);}
  }
  p->fn(p->data);
}
void GOMP_parallel(void (*fn)(void *),void *data,unsigned threads,unsigned flags) {
  parallel_fn real=(parallel_fn)dlsym(RTLD_NEXT,"GOMP_parallel");
  if(!real) _exit(127);
  Dl_info info;
  const char *caller=dladdr(__builtin_return_address(0),&info) ? info.dli_fname : "unknown";
  struct probe_data p={fn,data,threads,caller};
  real(probe_worker,&p,threads,flags);
}
