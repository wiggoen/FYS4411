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


double VariationalMonteCarlo::RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength)
{
    // Initialize matrices and variables
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double deltaEnergy;
    //double acceptanceWeight = 0;

    // Initial trial positions
    InitialTrialPositions(rOld, nParticles, nDimensions, stepLength);
    rNew = rOld;

    // Initialize classes
    Wavefunction waveFunction;
    Hamiltonian hamiltonian;


    // TODO: make own function for mcc?

    // loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        // Store the current value of the wave function
        waveFunctionOld = waveFunction.TrialWaveFunction(rOld, nParticles, nDimensions, alpha);
        waveFunction.QuantumForce(rOld, QForceOld, alpha);

        // New position to test
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < nDimensions; j++)
            {
                rNew(i, j) = rOld(i, j) + (RandomNumber() - 0.5) * stepLength;
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
            waveFunctionNew = waveFunction.TrialWaveFunction(rNew, nParticles, nDimensions, alpha);
            waveFunction.QuantumForce(rNew, QForceNew, alpha);

            // Sampling: Metropolis brute force
            MetropolisBruteForce(rNew, rOld, QForceOld, QForceNew, waveFunctionOld, waveFunctionNew, nDimensions, i);


            // update energies
            deltaEnergy = hamiltonian.LocalEnergy(rNew, nParticles, nDimensions, alpha);
            energySum += deltaEnergy;
            energySquaredSum += deltaEnergy*deltaEnergy;
        }
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    std::cout << "Energy: " << energy << "   &   Energy squared: " << energySquared << std::endl;
    return energy;
}


double VariationalMonteCarlo::RandomNumber()
{
    static std::random_device rd;  // Initialize the seed for the random number engine
    static std::mt19937_64 gen(rd());  // Call the Mersenne Twister algorithm
    // Set up the uniform distribution for x in [0, 1]
    static std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);
    // Set up the normal distribution for x in [0, 1]
    //std::normal_distribution<double> NormalDistribution(0.0,1.0);     // Will be used later
    return UniformNumberGenerator(gen);
}


void VariationalMonteCarlo::InitialTrialPositions(arma::mat &r, int nParticles, int nDimensions, int stepLength)
{
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            r(i, j) = (RandomNumber() - 0.5) * stepLength;
        }
    }
}


void VariationalMonteCarlo::MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew, int nDimensions, int i)
{
    // test is performed by moving one particle at the time
    // accept or reject this move
    double acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);

    //std::cout << acceptanceWeight << std::endl;

    if (RandomNumber() <= acceptanceWeight)
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
}
