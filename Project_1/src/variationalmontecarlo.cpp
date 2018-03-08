#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include <random>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h> /* Exit failure <- to force the program to stop: exit(EXIT_FAILURE);  */

VariationalMonteCarlo::VariationalMonteCarlo()
{

}


VariationalMonteCarlo::~VariationalMonteCarlo()
{

}


double VariationalMonteCarlo::RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength, int cycleStepToFile)
{
    // Adding variables to member variables
    this->nParticles = nParticles;
    this->nDimensions = nDimensions;
    this->nCycles = nCycles;
    this->alpha = alpha;
    this->stepLength = stepLength;
    this->cycleStepToFile = cycleStepToFile;

    // Initialize matrices and variables
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    x = 0;
    y = 0;
    acceptanceRatio = 0;

    // Initial trial positions
    InitialTrialPositions(rOld);
    rNew = rOld;

    // Run Monte Carlo cycles
    MonteCarloCycles();

    // Calculate energy
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    std::cout << std::setprecision(15) << "Energy: " << std::setw(5) << energy << "   &   Energy squared: " << energySquared << std::endl;
    std::cout << std::setprecision(15) << "Acceptance ratio = " << acceptanceRatio/(nCycles*nParticles)*100 << " %" << std::endl;
    std::cout << std::setprecision(15) << "Variance = " << (energySquared - energy*energy)/(nCycles * nParticles) << std::endl;
    return energy;
}


double VariationalMonteCarlo::RandomNumber()
{
    static std::random_device rd;  // Initialize the seed for the random number engine
    static std::mt19937_64 gen(rd());  // Call the Mersenne Twister algorithm
    // Set up the uniform distribution for x in [0, 1]
    static std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);
    return UniformNumberGenerator(gen);
}


double VariationalMonteCarlo::GaussianRandomNumber()
{
    static std::random_device rd;  // Initialize the seed for the random number engine
    static std::mt19937_64 gen(rd());  // Call the Mersenne Twister algorithm
    // Set up the normal distribution for x in [0, 1]
    static std::normal_distribution<double> NormalDistribution(0.0,1.0);
    return NormalDistribution(gen);
}


void VariationalMonteCarlo::InitialTrialPositions(arma::mat &r)
{
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            r(i, j) = (RandomNumber() - 0.5) * stepLength;
        }
    }
}


void VariationalMonteCarlo::MonteCarloCycles()
{
    // Initialize classes
    Wavefunction waveFunction;
    Hamiltonian hamiltonian;

    // Initialize variables
    waveFunctionOld = 0;
    waveFunctionNew = 0;
    energySum = 0;
    energySquaredSum = 0;
    deltaEnergy = 0;

    std::ofstream myfile;
    myfile.open("../Project_1/results.txt");

    // Loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        // Store the current value of the wave function
        waveFunctionOld = waveFunction.TrialWaveFunction(rOld, nParticles, nDimensions, alpha);
        waveFunction.QuantumForce(rOld, QForceOld, alpha);

        // New position to test
        for (int i = 0; i < nParticles; i++)
        {
            x = i; // Making i visible as private member x

            for (int j = 0; j < nDimensions; j++)
            {
                y = j; // Making j visible as private member y

                // Update position value
                rNew(x, y) = rOld(x, y) + (RandomNumber() - 0.5) * stepLength; // rOld(x, y) + (RandomNumber() - 0.5) * stepLength
                //std::cout << rNew << std::endl;
            }

            // Recalculate the value of the wave function and the quantum force
            waveFunctionNew = waveFunction.TrialWaveFunction(rNew, nParticles, nDimensions, alpha);
            waveFunction.QuantumForce(rNew, QForceNew, alpha);

            // Sampling: Metropolis brute force
            //MetropolisBruteForce(rNew, rOld, QForceOld, QForceNew, waveFunctionOld, waveFunctionNew);

            // Sampling: Fokker-Planck and Langevin
            FokkerPlanckAndLangevin(rNew, rOld, QForceOld, QForceNew, waveFunctionOld, waveFunctionNew);

            // Update energies
            deltaEnergy = hamiltonian.LocalEnergy(rNew, nParticles, nDimensions, alpha);
            energySum += deltaEnergy;
            energySquaredSum += deltaEnergy*deltaEnergy;
        }
        // Write to file
        if (cycle % cycleStepToFile == 0)
        {
            myfile << std::setw(10) << cycle << "     " << std::setprecision(6) << deltaEnergy << std::endl;
        }
    }
    myfile.close();
}


void VariationalMonteCarlo::MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew)
{
    double acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);

    // Test is performed by moving one particle at the time. Accept or reject this move.
    if (RandomNumber() <= acceptanceWeight)
    {
        rOld(x, y) = rNew(x, y);
        QForceOld(x, y) = QForceNew(x, y);
        waveFunctionOld = waveFunctionNew;
        acceptanceRatio += 1;
    } else
    {
        rNew(x, y) = rOld(x, y);
        QForceNew(x, y) = QForceOld(x, y);
    }

}

// TODO: IS THIS CORRECT IMPLEMENTED?
// acceptanceWeight around 1. Is this ok?
// Should we have used more abstraction in making the program? More like the mathematical equations?
void VariationalMonteCarlo::FokkerPlanckAndLangevin(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew)
{
    double D = 0.5;
    double dt = 0.01; // Interval [0.001,0.01]
    double acceptanceFactor = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
    double acceptanceWeight = (GreensFunction(rOld(x, y), rNew(x, y), D, dt, QForceOld(x, y))/GreensFunction(rNew(x, y), rOld(x, y), D, dt, QForceOld(x, y))) * acceptanceFactor;

    // Test is performed by moving one particle at the time. Accept or reject this move.
    if (RandomNumber() <= acceptanceWeight)
    {
        rOld(x, y) = rNew(x, y);
        QForceOld(x, y) = QForceNew(x, y);
        waveFunctionOld = waveFunctionNew;
        acceptanceRatio += 1;
    } else
    {
        rNew(x, y) = rOld(x, y) + D*QForceOld(x, y)*dt + GaussianRandomNumber()*sqrt(dt);
        QForceNew(x, y) = QForceOld(x, y);
    }
}

inline double VariationalMonteCarlo::GreensFunction(double oldPosition, double newPosition, double D, double dt, double QForceOld)
{
    return (1.0/pow(4.0*M_PI*D*dt, 3*nParticles/2.0)) * exp(-(newPosition-oldPosition-D*dt*QForceOld)*(newPosition-oldPosition-D*dt*QForceOld)/(4.0*D*dt)) + (nParticles - 1);
}
