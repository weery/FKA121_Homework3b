#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../rng_gen.h"
#include "../helper.h"

// -------------------------------
// DEFINES
// -------------------------------

// -------------------------------
// TASK SPECIFIC HELPER FUNCTIONS
// -------------------------------

// -------------------------------
// MAIN PROGRAM
// -------------------------------
int main()
{
    // -------------------------------
    // Initialize the gsl random number generator
    // -------------------------------
    Initialize_Generator();

    // -------------------------------
    // Variable Declarations
    // -------------------------------
    FILE* file;

    // -------------------------------
    // Initialize Variables
    // -------------------------------

    // -------------------------------
    // Simulation
    // -------------------------------


    // -------------------------------
    // Print results to file
    // -------------------------------
    file = fopen("*.dat","w");
    for (int i = 0; i < 10; ++i)
    {
        fprintf(file, "%i\n", i);
    }
    fclose(file);

    // -------------------------------
    // Free the gsl random number generator
    // -------------------------------
    Free_Generator();

    return 0;
}

// -------------------------------
// HELPER FUNCTION DEFINITIONS
// -------------------------------
