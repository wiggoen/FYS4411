#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include <random>
#include <iomanip>
#include <iostream>
#include <armadillo>

VariationalMonteCarlo::VariationalMonteCarlo() : // TODO: put as command line arguments
    nDimensions(1),
    nParticles(1),
    nCycles(1e6),
    alpha(0.5),
    stepLength(0.1)
{

}

VariationalMonteCarlo::~VariationalMonteCarlo()
{

}


void VariationalMonteCarlo::runMonteCarloIntegration()
{
    // TODO: Move random to main?

    // Initialize the seed and call the Mersenne Twister algorithm
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x in [0, 1]
    std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);
    // Set up the normal distribution for x in [0, 1]
    //std::normal_distribution<double> NormalDistribution(0.0,1.0);         // Will be used later


    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaEnergy;

    double acceptanceWeight = 0;

    // initial trial positions
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rOld(i, j) = (UniformNumberGenerator(gen) - 0.5) * stepLength;
        }
    }
    rNew = rOld;

    // TODO: make own function for mcc?

    // loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        // Store the current value of the wave function
        waveFunctionOld = waveFunction(rOld);
        QuantumForce(rOld, QForceOld);

        // New position to test
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < nDimensions; j++)
            {
                rNew(i, j) = rOld(i, j) + (UniformNumberGenerator(gen) - 0.5) * stepLength;
            }
            //  for the other particles we need to set the position to the old position since
            //  we move only one particle at the time
            for (int k = 0; k < nParticles; k++) {
                if ( k != i)
                {
                    for (int j = 0; j < nDimensions; j++)
                    {
                        rNew(k, j) = rOld(k, j);
                    }
                }
            }
            // Recalculate the value of the wave function and the quantum force
            waveFunctionNew = waveFunction(rNew);
            QuantumForce(rNew, QForceNew);

            // TODO: Move sampling

            // Metropolis brute force
            // test is performed by moving one particle at the time
            // accept or reject this move
            acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);

            //std::cout << acceptanceWeight << std::endl;

            if (UniformNumberGenerator(gen) <= acceptanceWeight)
            {
                for (int j = 0; j < nDimensions; j++)
                {
                    rOld(i, j) = rNew(i, j);
                    QForceOld(i, j) = QForceNew(i, j);
                    waveFunctionOld = waveFunctionNew;
                }
            } else
            {
                for (int j = 0; j < nDimensions; j++)
                {
                    rNew(i, j) = rOld(i, j);
                    QForceNew(i, j) = QForceOld(i, j);
                }
            }
            // update energies
            deltaEnergy = localEnergy(rNew);
            energySum += deltaEnergy;
            energySquaredSum += deltaEnergy*deltaEnergy;
        }
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    std::cout << "Energy: " << energy << "   &   Energy squared: " << energySquared << std::endl;
}
