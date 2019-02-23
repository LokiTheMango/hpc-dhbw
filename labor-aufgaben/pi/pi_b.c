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
  int globalCount = 0, globalSamples = TRYS;

#pragma omp parallel shared(globalSamples, globalCount)
  {
    int lul = 0;
#pragma omp for
    for (int i = 0; i < globalSamples; ++i)
    {
      lul += throw();
    }

#pragma omp critical
    globalCount += lul;
  }

  double pi = 4.0 * (double)globalCount / (double)(globalSamples);

  printf("pi is %.9lf\n", pi);

  return 0;
}
