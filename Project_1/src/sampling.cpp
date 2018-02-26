#include "inc/sampling.h"


Sampling::Sampling()
{

}

Sampling::~Sampling()
{

}

// TODO: FIX THIS
/*
Sampling::Metropolis_brute_force(double *matrix)
{
    // Loop over matrix
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            // Find random position
            double x = (double) RandomNumberGenerator(gen);
            double y = (double) RandomNumberGenerator(gen);

            // Compute deltaE
            double deltaEnergy = computeDeltaE(x, y);


            // Metropolis test
            if (RandomNumberGenerator(gen) <= SOME_VALUE)
            {
                matrix[x][y] = 0;                // Change something
                totalEnergy += deltaEnergy;         // Update energy
            }
        }
    }
    return matrix;
}
*/
