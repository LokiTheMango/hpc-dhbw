#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  int nthreads, tid;

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel private(nthreads, tid)
  {

    /* Obtain thread number */
    nthreads = omp_get_num_threads();
    tid = omp_get_thread_num();
    printf("Hello World from thread = %d of %d\n", tid, nthreads);
  }
}
