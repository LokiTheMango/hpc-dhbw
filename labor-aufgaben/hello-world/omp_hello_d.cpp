#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  int nthreads, tid;

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel sections
  {
    printf("Hallo Welt\n");
#pragma omp section
    printf("Hello World\n");
#pragma omp section
    printf("Bonjour monde\n");
#pragma omp section
    printf("OOOF monde\n");
  }
}
