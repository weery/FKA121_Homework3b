#ifndef _helper_h_
#define _helper_h_

// -------------------------------
// Includes required for header to compile
// -------------------------------
#include <complex.h>

// -------------------------------
// Defines
// -------------------------------
#define HBAR    0.1
#define PI      3.14159

// -------------------------------
// Function declarations
// -------------------------------
double complex Gaussian_Wave_Packet (double, double, double, double);

void Wave_To_Probability_Density (double*, double*, unsigned int);

#endif
