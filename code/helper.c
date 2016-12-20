#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "./helper.h"
#include <complex.h>

// -------------------------------
// Function definitions
// -------------------------------
double complex Gaussian_Wave_Packet (double x, double d, double x_0, double p_0)
{
    double prefactor            = pow(1/(PI*d*d),1.0/4);
    double exp_1                = -(x-x_0)*(x-x_0)/(2*d*d);
    double complex exp_2        = I * p_0 * (x-x_0)/HBAR;

    return prefactor*cexp(exp_1+exp_2);
}

void Wave_To_Probability_Density (double* denisty, double* wave, unsigned int n)
{
    for (int i = 0; i < n; i++)
    {
        density[i] = wave[i]*wave[i];
    }
}

void Position_To_Moment_Space(double * moment, double* position, unsigned int n)
{
    
}
