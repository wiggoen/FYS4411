#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include <random>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>


VariationalMonteCarlo::VariationalMonteCarlo() {}


VariationalMonteCarlo::~VariationalMonteCarlo() {}


double VariationalMonteCarlo::RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength)
{
    // Adding variables to member variables
    this->nParticles = nParticles;
    this->nDimensions = nDimensions;
    this->nCycles = nCycles;
    this->alpha = alpha;
    this->stepLength = stepLength;

    // Initialize matrices and variables
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);

    // Initial trial positions
    InitialTrialPositions(rOld);
    rNew = rOld;

    // Run Monte Carlo cycles
    MonteCarloCycles();

    // Calculate energy
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
            for (int j = 0; j < nDimensions; j++)

            {
                rNew(i, j) = rOld(i, j) + (RandomNumber() - 0.5) * stepLength;
            }
            // For the other particles we need to set the position to the old position since we move only one particle at the time
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
            //MetropolisBruteForce(rNew, rOld, QForceOld, QForceNew, waveFunctionOld, waveFunctionNew, i);

            // Sampling: Fokker-Planck and Langevin
            FokkerPlanckAndLangevin(rNew, rOld, QForceOld, QForceNew, waveFunctionOld, waveFunctionNew, i);

            // Update energies
            deltaEnergy = hamiltonian.LocalEnergy(rNew, nParticles, nDimensions, alpha);
            energySum += deltaEnergy;
            energySquaredSum += deltaEnergy*deltaEnergy;

        }
        // Write to file
        int num = 1000;
        if (cycle % num == 0)
        {
            myfile << std::setw(10) << cycle << "     " << std::setprecision(6) << deltaEnergy << std::endl;
        }
    }
    myfile.close();
}


void VariationalMonteCarlo::MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew, int i)
{
    double acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);

    // Test is performed by moving one particle at the time. Accept or reject this move.
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

// TODO: Fix FokkerPlank importance sampling
// acceptanceWeight around 1. Is this ok?
// Should we have used more abstraction in making the program? More like the mathematical equations?
void VariationalMonteCarlo::FokkerPlanckAndLangevin(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew, int i)
{
    double acceptanceFactor = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
    double D = 0.5;
    double deltaT = 0.01; // Interval [0.001,0.01]
    double acceptanceWeight = 0;

    // Need if-test and changes in matrices for rNew, rOld, QForceOld, QForceNew
    for (int j = 0; j < nDimensions; j++)
    {
        acceptanceWeight = (GreensFunction(rOld(i, j), rNew(i, j), D, deltaT, QForceOld(i, j))/GreensFunction(rNew(i, j), rOld(i, j), D, deltaT, QForceOld(i, j))) * acceptanceFactor;
        //std::cout << "acceptanceWeight = " << acceptanceWeight << std::endl;
    }

    // Test is performed by moving one particle at the time. Accept or reject this move.
    //if (RandomNumber() <= acceptanceWeight)
    /*
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
            rNew(i, j) = rOld(i, j) + D*QForceOld(i, j)*deltaT + GaussianRandomNumber()*sqrt(deltaT);
            QForceNew(i, j) = QForceOld(i, j);
        }
    }
    */
}

// TODO: Fix nParticles dependence. Should take nParticles from class.

double VariationalMonteCarlo::GreensFunction(double x, double y, double D, double deltaT, double QForceOld)
{
    return (1.0/pow(4.0*M_PI*D*deltaT, 3*nParticles/2)) * exp(-(y-x-D*deltaT*QForceOld)*(y-x-D*deltaT*QForceOld)/(4.0*D*deltaT));
}
