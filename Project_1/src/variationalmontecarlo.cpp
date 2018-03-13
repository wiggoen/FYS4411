#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>  // high resolution timing: http://en.cppreference.com/w/cpp/chrono/c/clock
#include <math.h>
//#include <stdlib.h> /* Exit failure <- to force the program to stop: exit(EXIT_FAILURE);  */


VariationalMonteCarlo::VariationalMonteCarlo()
{

}


VariationalMonteCarlo::~VariationalMonteCarlo()
{

}


arma::rowvec VariationalMonteCarlo::RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles,
                                                             double alpha, double stepLength, double dt,
                                                             int cycleStepToFile)
{
    // Adding variables to member variables
    this->nParticles = nParticles;
    this->nDimensions = nDimensions;
    this->nCycles = nCycles;
    this->alpha = alpha;
    this->stepLength = stepLength;
    this->dt = dt;
    this->cycleStepToFile = cycleStepToFile;

    // Initialize matrices and variables
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    energySum = 0;
    energySquaredSum = 0;
    deltaEnergy = 0;
    acceptanceCounter = 0;

    arma::rowvec runDetails;
    double time = 0;
    double energy = 0;
    double energySquared = 0;
    double variance = 0;
    double acceptanceRatio = 0;

    // Initial trial positions
    //InitialTrialPositionsBruteForce(rOld);
    InitialTrialPositionsImportanceSampling(rOld);
    rNew = rOld;

    // Store the current value of the wave function and quantum force
    waveFunctionOld = Wavefunction::TrialWaveFunction(rOld, nParticles, nDimensions, alpha);
    Wavefunction::QuantumForce(rOld, QForceOld, alpha);
    QForceNew = QForceOld;

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Run Monte Carlo cycles
    MonteCarloCycles();

    // Timing finished
    auto end_time = std::chrono::high_resolution_clock::now();
    time = (std::chrono::duration<double> (end_time - start_time).count());

    // Calculation
    energy = energySum/(nCycles * nParticles);
    energySquared = energySquaredSum/(nCycles * nParticles);

    variance = (energySquared - energy*energy)/(nCycles * nParticles);
    acceptanceRatio = acceptanceCounter/(nCycles * nParticles);

    runDetails << time << energy << energySquared << variance << acceptanceRatio;
    return runDetails;
}


void VariationalMonteCarlo::InitialTrialPositionsBruteForce(arma::mat &r)
{
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            r(i, j) = (UniformRandomNumber() - 0.5) * stepLength;
        }
    }
}


void VariationalMonteCarlo::InitialTrialPositionsImportanceSampling(arma::mat &r)
{
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            r(i, j) = GaussianRandomNumber()*sqrt(dt);
        }
    }
}


void VariationalMonteCarlo::MonteCarloCycles()
{
    // Initialize variables
    waveFunctionOld = 0;
    waveFunctionNew = 0;

    std::ofstream myfile;
    myfile.open("../Project_1/results.txt");

    // Loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        // Sampling
        //MetropolisBruteForce(rNew, rOld, waveFunctionOld, waveFunctionNew);
        ImportanceSampling(rNew, rOld, QForceOld, QForceNew, waveFunctionOld, waveFunctionNew);

        // Write to file
        if (cycle % cycleStepToFile == 0)
        {
            myfile << std::setw(10) << cycle << "     " << std::setprecision(6) << deltaEnergy << std::endl;
        }
    }
    myfile.close();
}


double VariationalMonteCarlo::UniformRandomNumber()
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


void VariationalMonteCarlo::MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, double &waveFunctionOld,
                                                 double &waveFunctionNew)
{
    // New position to test
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rNew(i, j) = rOld(i, j) + (UniformRandomNumber() - 0.5) * stepLength;
        }

        // Recalculate the value of the wave function
        waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, nParticles, nDimensions, alpha);

        acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);

        UpdateEnergies(i);
    }
}


void VariationalMonteCarlo::ImportanceSampling(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld,
                                               arma::mat &QForceNew, double &waveFunctionOld,
                                               double &waveFunctionNew)
{
    double D = 0.5; // Diffusion coefficient
    double acceptanceFactor = 0;
    double GreensOldNew = 0;
    double GreensNewOld = 0;

    // New position to test
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rNew(i, j) = rOld(i, j) + D*QForceOld(i, j)*dt + GaussianRandomNumber()*sqrt(dt);
        }

        // Recalculate the value of the wave function and the quantum force
        waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, nParticles, nDimensions, alpha);
        Wavefunction::QuantumForce(rNew, QForceNew, alpha);

        acceptanceFactor = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
        GreensOldNew = GreensFunction(rOld, rNew, QForceOld, D, dt, i);
        GreensNewOld = GreensFunction(rNew, rOld, QForceNew, D, dt, i);
        acceptanceWeight = (GreensOldNew/GreensNewOld) * acceptanceFactor;

        UpdateEnergies(i);
    }
}


double VariationalMonteCarlo::GreensFunction(const arma::mat &rOld, const arma::mat &rNew, const arma::mat &QForceOld,
                                             double &D, double &dt, int &i)
{
    double fourDdt = 4.0*D*dt;
    arma::rowvec yx = rNew.row(i)-rOld.row(i)-D*dt*QForceOld.row(i);
    double yxSquared = arma::dot(yx, yx);
    return exp(-yxSquared/fourDdt) + (nParticles - 1);
}


void VariationalMonteCarlo::UpdateEnergies(int &i)
{
    // Test is performed by moving one particle at the time. Accept or reject this move.
    if (UniformRandomNumber() <= acceptanceWeight)
    {
        rOld.row(i) = rNew.row(i);
        QForceOld.row(i) = QForceNew.row(i);
        waveFunctionOld = waveFunctionNew;
        acceptanceCounter += 1;
    } else
    {
        rNew.row(i) = rOld.row(i);
        QForceNew.row(i) = QForceOld.row(i);
    }

    // Update energies
    deltaEnergy = Hamiltonian::LocalEnergy(rNew, nParticles, nDimensions, alpha);
    energySum += deltaEnergy;
    energySquaredSum += deltaEnergy*deltaEnergy;
}

void VariationalMonteCarlo::Blocking(arma::vec samples, int nSamples, int block_size, arma::vec results)
{
    //Integer division will waste some samples
    int nBlocks = nSamples/block_size;
    double block_samples = new arma::vec(nBlocks);
}

double VariationalMonteCarlo::Mean(arma::vec &samples, int nSamples)
{
    double m = 0;
    for(int i=0; i<nSamples; i++)
    {
        m+=samples(i);
    }
    return m/double(nSamples);
}

double VariationalMonteCarlo::Variance(arma::vec &samples, int nSamples, arma::vec &results)
{
    double m2=0, m=0;
    for(int i=0; i<nSamples; i++)
    {
        m+=samples(i);
        m2+=samples(i)*samples(i);
    }
    m /= double(nSamples);
    m2 /= double(nSamples);
    results(0) = m;
    results(1) = m2-(m*m);
}
