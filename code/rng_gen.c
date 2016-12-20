#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static gsl_rng *q;

void Initialize_Generator()
{
	gsl_rng_env_setup();
    const gsl_rng_type *T;
    T = gsl_rng_default;
    q = gsl_rng_alloc(T);
    gsl_rng_set(q,time(NULL));
}

double rand_uniform()
{
    return gsl_rng_uniform(q);
}

void Free_Generator()
{
    gsl_rng_free(q);
}

int randint(int min, int max)
{
	double r= rand_uniform();
	int val = (int)(min +(double)(max-min)*r+0.5);
	return val;
}

double randdouble(double min, double max)
{
	double r = rand_uniform();
	double val = min+(max-min)*r;
	return val;
}
