#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define TRIALS 1e8

int main() {
  float* x = malloc(TRIALS*sizeof(float));
  float* y = malloc(TRIALS*sizeof(float));
  int inside = 0;

  /* Seed RNG with current time */
  srand((unsigned)time(NULL)); 

  /* Generate random points inside a centered square of area 4
     and count the proportion that falls within the unit circle. */
  int i;
  for (i = 0; i < TRIALS; i++) {
    x[i] = 2.0 * (rand()/(RAND_MAX+1.0)) - 1.0;
    y[i] = 2.0 * (rand()/(RAND_MAX+1.0)) - 1.0;
    if (sqrt(x[i]*x[i]+y[i]*y[i]) < 1.0) {
      inside++;
    }
  }
 
  free(x);
  free(y);

  printf("Pi is approximately %.6f.\n", 4.0*inside/TRIALS); 

  return 0;
}
