#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>    /* high resolution timing: http://en.cppreference.com/w/cpp/chrono/c/clock */
#include <cmath>
#include <stdlib.h>  /* Exit failure, to force the program to stop: exit(EXIT_FAILURE); */


VariationalMonteCarlo::VariationalMonteCarlo() :
    nDimensions(2)
{

}


VariationalMonteCarlo::~VariationalMonteCarlo( void )
{

}

arma::rowvec VariationalMonteCarlo::RunMonteCarloIntegration(const int nParticles, const int nCycles,
                                                             const double alpha, const double beta, const double omega, const double a, const double stepLength)
{
    /* Adding variables to member variables */
    this->nParticles      = nParticles;
    this->nCycles         = nCycles;
    this->alpha           = alpha;
    this->beta            = beta;
    this->omega           = omega;
    this->stepLength      = stepLength;
    //this->timeStep        = timeStep;
    //this->cycleStepToFile = cycleStepToFile;
    //this->samplingType    = samplingType;
    //this->derivationType  = derivationType;
    //this->cycleType       = cycleType;

    samplingType = "BruteForce";
    cycleType = "MonteCarlo";
    derivationType = "Analytical";

    /* Initialize matrices and variables */
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    /*
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions); */
    waveFunctionOld   = 0.0;
    waveFunctionNew   = 0.0;
    energySum         = 0.0;
    energySquaredSum  = 0.0;
    deltaEnergy       = 0.0;
    acceptanceWeight  = 0.0;
    acceptanceCounter = 0.0;
    /*
    psiSum            = 0.0;
    psiTimesEnergySum = 0.0;
    deltaPsi          = 0.0; */


    arma::rowvec runDetails;
    double runTime = 0.0;
    double energy = 0.0;
    double energySquared = 0.0;
    double variance = 0.0;
    double acceptanceRatio = 0.0;

    /*
    if (derivationType == "Interaction")  { beta = 2.82843; }
    else                                  { beta = 1.0; }
    */
    /* Initial trial positions */
    if (samplingType == "BruteForce")      { InitialTrialPositionsBruteForce(rOld); }
    //else if (samplingType == "Importance") { InitialTrialPositionsImportanceSampling(rOld); }

    //if (derivationType == "Interaction")   { CheckInitialDistance(rOld); }

    rNew = rOld;

    /* Store the current value of the wave function and quantum force */
    waveFunctionOld = Wavefunction::TrialWaveFunction(rOld, alpha, beta, omega, a, 1);
    //Wavefunction::QuantumForce(rOld, QForceOld, alpha);
    //QForceNew = QForceOld;

    /* Start timing */
    auto start_time = std::chrono::high_resolution_clock::now();

    /* Run Monte Carlo cycles */
    if (cycleType == "MonteCarlo")//|| cycleType == "OneBodyDensity")
    {
        std::cout << "Running MC Cycles.." << std::endl;
        MonteCarloCycles();
    } /*else if (cycleType == "SteepestDescent") {
        SteepestDescent(nParticles);
    }*/

    /* Timing finished */
    auto end_time = std::chrono::high_resolution_clock::now();
    runTime = (std::chrono::duration<double> (end_time - start_time).count());

    /* Normalizing */
    double normalizationFactor = 1.0/(nCycles * nParticles);

    /* Calculation of averages */
    energy = energySum * normalizationFactor;
    energySquared = energySquaredSum * normalizationFactor;

    variance = (energySquared - energy*energy);
    acceptanceRatio = acceptanceCounter * normalizationFactor;

    /*
    if (cycleType == "OneBodyDensity")
    {
        for (int i = 0; i < nBins; i++)
        {
            hist(i) /= volume(i) * normalizationFactor;
        }
        std::ofstream histogram;  // Write position matrix to file
        histogram.open("../Project_1/results/histogram.txt");

        // Loop over Monte Carlo cycles
        for (int i = 0; i < nBins; i++)
        {
            // Write to file
            histogram << hist(i) << std::endl;
        }
        histogram.close();
    }
    */

    runDetails << runTime << energy << energySquared << variance << acceptanceRatio;
    return runDetails;
}


void VariationalMonteCarlo::InitialTrialPositionsBruteForce(arma::mat &r)
{
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            r(i, j) = (UniformRandomNumber() - 0.5);
        }
    }
}


void VariationalMonteCarlo::MonteCarloCycles( void )
{
    //std::ofstream myfile;
    //myfile.open("../Project_1/results/results.txt");

    /* Loop over Monte Carlo cycles */
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        /* Sampling */
        if (samplingType == "BruteForce") {
            MetropolisBruteForce(rNew, rOld, waveFunctionNew, waveFunctionOld);
        } /*else if (samplingType == "Importance")
        {
            ImportanceSampling(rNew, rOld, QForceNew, QForceOld, waveFunctionNew, waveFunctionOld);
        } */

        /* Write to file */
        /*
        if (cycle % cycleStepToFile == 0)
        {
            myfile << std::setw(10) << cycle << "     " << std::setprecision(6) << energySum/(cycle*nParticles) << std::endl;
        } */
    }
    //myfile.close();
}


double VariationalMonteCarlo::UniformRandomNumber( void )
{
    static std::random_device rd;  /* Initialize the seed for the random number engine */
    static std::mt19937_64 gen(rd());  /* Call the Mersenne Twister algorithm */
    /* Set up the uniform distribution for x in [0, 1] */
    static std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);
    return UniformNumberGenerator(gen);
}


void VariationalMonteCarlo::MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, double &waveFunctionNew,
                                                 const double &waveFunctionOld)
{
    /* New position to test */
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rNew(i, j) = rOld(i, j) + (UniformRandomNumber() - 0.5) * stepLength;
        }

        /* Recalculate the value of the wave function */
        /*if (derivationType == "Interaction")
        {
            waveFunctionNew = Wavefunction::TrialWaveFunctionInteraction(rNew, nParticles, nDimensions, alpha, beta, a);
        } else
        {*/
        waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, alpha, beta, omega, a, 1);
        //}

        acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
        //std::cout << "wOld = " << waveFunctionOld << std::endl;
        //std::cout << "wNew = " <<waveFunctionNew << std::endl;
        //std::cout << "acceptW = " <<acceptanceWeight << std::endl;

        UpdateEnergies(i);
    }
}


void VariationalMonteCarlo::UpdateEnergies(const int &i)
{
    /* Test is performed by moving one particle at the time. Accept or reject this move. */
    if (UniformRandomNumber() <= acceptanceWeight)
    {
        rOld.row(i) = rNew.row(i);
        //QForceOld.row(i) = QForceNew.row(i);
        waveFunctionOld = waveFunctionNew;
        acceptanceCounter += 1;
    } else
    {
        rNew.row(i) = rOld.row(i);
        //QForceNew.row(i) = QForceOld.row(i);
        //waveFunctionNew = waveFunctionOld;  /* Probably unnecessary since the wavefunction didn't change. */
    }

    if (derivationType == "Analytical") {
        /* Update energies (without numerical derivation and interaction) */
        deltaEnergy = Hamiltonian::LocalEnergy(rNew, nParticles, alpha, beta, omega, a);
    } /*else if (derivationType == "Numerical")
    {
        // Update energies using numerical derivation
        deltaEnergy = Hamiltonian::NumericalLocalEnergy(rNew, nParticles, nDimensions, alpha, stepLength, beta);
    } else if (derivationType == "Interaction")
    {
        // Update energies using interaction
        deltaEnergy = Hamiltonian::LocalEnergyInteraction(rNew, nParticles, nDimensions, alpha, beta, a);
    }*/

    energySum         += deltaEnergy;
    energySquaredSum  += (deltaEnergy*deltaEnergy);

    /*
    if (cycleType == "OneBodyDensity")
    {
        OneBodyDensity();
    }


    if (cycleType == "SteepestDescent")
    {
        deltaPsi           = Wavefunction::DerivativePsi(rNew, nParticles, nDimensions, beta);
        psiSum            += deltaPsi;
        psiTimesEnergySum += deltaPsi*deltaEnergy;
    }*/
}
