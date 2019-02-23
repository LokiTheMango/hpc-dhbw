#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>

#define TRYS 5000000

static int throw()
{
  double x, y;
  x = (double)rand() / (double)RAND_MAX;
  y = (double)rand() / (double)RAND_MAX;
  if ((x * x + y * y) <= 1.0)
    return 1;

  return 0;
}

int main(int argc, char **argv)
{
  // Steuern: env -> OMP_NUM_THREADS=6
  // Unterbinden
  omp_set_num_threads(6);

  int globalCount = 0, globalSamples = TRYS;

#pragma omp parallel shared(globalSamples) reduction(+ \
                                                     : globalCount)
  {
#pragma omp for
    for (int i = 0; i < globalSamples; ++i)
    {
      globalCount += throw();
    }

#pragma omp critical
    printf("Thread %d treffer: %d\n", omp_get_thread_num(), globalCount);
  }

  double pi = 4.0 * (double)globalCount / (double)(globalSamples);

  printf("pi is %.9lf\n", pi);

  return 0;
}
