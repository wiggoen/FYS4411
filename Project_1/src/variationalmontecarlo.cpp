#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include <random>
#include <iomanip>
#include <iostream>


VariationalMonteCarlo::VariationalMonteCarlo()
{

}

VariationalMonteCarlo::~VariationalMonteCarlo()
{

}


double VariationalMonteCarlo::runMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength)
{
    std::random_device rd;                                                   // Initialize the seed
    std::mt19937_64 gen(rd());                                               // Call the Mersenne Twister algorithm
    std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);  // Set up the uniform distribution for x in [0, 1]
    //std::normal_distribution<double> NormalDistribution(0.0,1.0);            // Set up the normal distribution for x in [0, 1]


    // Initialize matrices and variables
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);

    double waveFunctionOld = 0; // TODO: move to wavefunction?
    double waveFunctionNew = 0;

    double energySum = 0; // TODO: move to hamiltonian?
    double energySquaredSum = 0;

    double deltaEnergy;

    double acceptanceWeight = 0;

    Wavefunction waveFunction;
    Hamiltonian hamiltonian;

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
        waveFunctionOld = waveFunction.trialWaveFunction(rOld, nParticles, nDimensions, alpha);
        waveFunction.QuantumForce(rOld, QForceOld, alpha);

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
            waveFunctionNew = waveFunction.trialWaveFunction(rNew, nParticles, nDimensions, alpha);
            waveFunction.QuantumForce(rNew, QForceNew, alpha);



            // TODO: Move sampling to own function

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
            deltaEnergy = hamiltonian.localEnergy(rNew, nParticles, nDimensions, alpha);
            energySum += deltaEnergy;
            energySquaredSum += deltaEnergy*deltaEnergy;
        }
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    std::cout << "Energy: " << energy << "   &   Energy squared: " << energySquared << std::endl;
    return energy;
}
