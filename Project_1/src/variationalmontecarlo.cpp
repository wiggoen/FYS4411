#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>    // high resolution timing: http://en.cppreference.com/w/cpp/chrono/c/clock
#include <cmath>
#include <stdlib.h>  // Exit failure, to force the program to stop: exit(EXIT_FAILURE);


VariationalMonteCarlo::VariationalMonteCarlo() :
    a(0.0043)  // hard-core diameter of the bosons
{

}


VariationalMonteCarlo::~VariationalMonteCarlo()
{

}


arma::rowvec VariationalMonteCarlo::RunMonteCarloIntegration(const int nParticles, const int nDimensions, const int nCycles,
                                                             const double alpha, const double stepLength,
                                                             const double timeStep, const int cycleStepToFile,
                                                             const std::string samplingType,
                                                             const std::string derivationType,
                                                             const std::string cycleType)
{
    // Adding variables to member variables
    this->nParticles      = nParticles;
    this->nDimensions     = nDimensions;
    this->nCycles         = nCycles;
    this->alpha           = alpha;
    this->stepLength      = stepLength;
    this->timeStep        = timeStep;
    this->cycleStepToFile = cycleStepToFile;
    this->samplingType    = samplingType;
    this->derivationType  = derivationType;
    this->cycleType       = cycleType;

    // Initialize matrices and variables
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    waveFunctionOld   = 0.0;
    waveFunctionNew   = 0.0;
    energySum         = 0.0;
    energySquaredSum  = 0.0;
    deltaEnergy       = 0.0;
    acceptanceWeight  = 0.0;
    acceptanceCounter = 0.0;
    thrownCounter     = 0.0;
    psiSum            = 0.0;
    psiTimesEnergySum = 0.0;
    deltaPsi          = 0.0;

    arma::rowvec runDetails;
    double runTime = 0.0;
    double energy = 0.0;
    double energySquared = 0.0;
    double variance = 0.0;
    double acceptanceRatio = 0.0;

    if (derivationType == "Interaction")  { beta = 2.82843; }
    else                                  { beta = 1.0; }

    // Initial trial positions
    if (samplingType == "BruteForce")      { InitialTrialPositionsBruteForce(rOld); }
    else if (samplingType == "Importance") { InitialTrialPositionsImportanceSampling(rOld); }

    if (derivationType == "Interaction")   { CheckInitialDistance(rOld); }

    rNew = rOld;

    // Store the current value of the wave function and quantum force
    waveFunctionOld = Wavefunction::TrialWaveFunction(rOld, nParticles, nDimensions, alpha, beta);
    Wavefunction::QuantumForce(rOld, QForceOld, alpha);
    QForceNew = QForceOld;

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Run Monte Carlo cycles
    if (cycleType == "MonteCarlo") {
        MonteCarloCycles();
    } else if (cycleType == "SteepestDescent") {
        SteepestDescent(nParticles);
    }

    // Timing finished
    auto end_time = std::chrono::high_resolution_clock::now();
    runTime = (std::chrono::duration<double> (end_time - start_time).count());

    // Normalizing
    int nCyclesNew = nCycles - thrownCounter;
    double normalizationFactor = 1.0/(nCyclesNew * nParticles);

    // Calculation
    energy = energySum * normalizationFactor;
    energySquared = energySquaredSum * normalizationFactor;

    variance = (energySquared - energy*energy);
    acceptanceRatio = acceptanceCounter * normalizationFactor;

    runDetails << runTime << energy << energySquared << variance << acceptanceRatio;
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
            r(i, j) = GaussianRandomNumber()*sqrt(timeStep);
        }
    }
}


void VariationalMonteCarlo::RedrawPositionImportanceSampling(arma::mat &r, const int &i)
{
    for (int j = 0; j < nDimensions; j++)
    {
        if (derivationType == "Importance")
        {
            r(i, j) = GaussianRandomNumber()*sqrt(timeStep);
        } else
        {
            r(i, j) = (UniformRandomNumber() - 0.5) * stepLength;
        }
    }
}


void VariationalMonteCarlo::CheckInitialDistance(arma::mat &rOld)
{
    double distance = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            distance = Hamiltonian::ParticleDistance(rOld.row(i), rOld.row(j));
            while (distance < a)
            {
                // Making sure that no particle lies within a distance 'a' from each other
                RedrawPositionImportanceSampling(rOld, i);
                distance = Hamiltonian::ParticleDistance(rOld.row(i), rOld.row(j));
            }
        }
    }
}


void VariationalMonteCarlo::MonteCarloCycles()
{
    std::ofstream myfile;
    myfile.open("../Project_1/results.txt");

    // Loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        // Sampling
        if (samplingType == "BruteForce") {
            MetropolisBruteForce(rNew, rOld, waveFunctionNew, waveFunctionOld);
        } else if (samplingType == "Importance")
        {
            ImportanceSampling(rNew, rOld, QForceNew, QForceOld, waveFunctionNew, waveFunctionOld);
        }

        // Write to file
        if (cycle % cycleStepToFile == 0)
        {
            myfile << std::setw(10) << cycle << "     " << std::setprecision(6) << energySum/(cycle*nParticles) << std::endl;
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


void VariationalMonteCarlo::MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, double &waveFunctionNew,
                                                 const double &waveFunctionOld)
{
    // New position to test
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rNew(i, j) = rOld(i, j) + (UniformRandomNumber() - 0.5) * stepLength;
        }

        // Recalculate the value of the wave function
        if (derivationType == "Interaction")
        {
            waveFunctionNew = Wavefunction::TrialWaveFunctionInteraction(rNew, nParticles, nDimensions, alpha, beta, a);
        } else
        {
            waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, nParticles, nDimensions, alpha, beta);
        }

        acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);

        UpdateEnergies(i);
    }
}


void VariationalMonteCarlo::ImportanceSampling(arma::mat &rNew, const arma::mat &rOld, arma::mat &QForceNew,
                                               const arma::mat &QForceOld, double &waveFunctionNew,
                                               const double &waveFunctionOld)
{
    double diffusionCoefficient = 0.5;
    double wavefunctionsSquared = 0.0;
    double GreensRatio = 0.0;

    // New position to test
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rNew(i, j) = rOld(i, j) + diffusionCoefficient*QForceOld(i, j)*timeStep + GaussianRandomNumber()*sqrt(timeStep);
        }

        // Recalculate the value of the wave function and the quantum force
        if (derivationType == "Analytical")
        {
            waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, nParticles, nDimensions, alpha, beta);
            Wavefunction::QuantumForce(rNew, QForceNew, alpha);
        } else if (derivationType == "Numerical")
        {
            waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, nParticles, nDimensions, alpha, beta);
            Wavefunction::NumericalQuantumForce(rNew, QForceNew, nParticles, nDimensions, alpha, stepLength, beta);
        } else if (derivationType == "Interaction")
        {
            waveFunctionNew = Wavefunction::TrialWaveFunctionInteraction(rNew, nParticles, nDimensions, alpha, beta, a);
            Wavefunction::QuantumForceInteraction(rNew, QForceNew, nParticles, nDimensions, alpha, beta, a, i);
        }

        wavefunctionsSquared = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
        GreensRatio = GreensFunction(rNew, rOld, QForceNew, QForceOld, diffusionCoefficient, timeStep, i);
        acceptanceWeight = GreensRatio * wavefunctionsSquared;

        UpdateEnergies(i);
    }
}


double VariationalMonteCarlo::GreensFunction(const arma::mat &rNew, const arma::mat &rOld, const arma::mat &QForceNew,
                                             const arma::mat &QForceOld, const double &diffusionCoefficient,
                                             const double &timeStep, const int &i)
{
    double fraction = 1.0/(4.0*diffusionCoefficient*timeStep);
    arma::rowvec yx = rNew.row(i) - rOld.row(i) - diffusionCoefficient*timeStep*QForceOld.row(i);
    double yxSquared = arma::dot(yx, yx);
    arma::rowvec xy = rOld.row(i) - rNew.row(i) - diffusionCoefficient*timeStep*QForceNew.row(i);
    double xySquared = arma::dot(xy, xy);
    return exp((-xySquared + yxSquared) * fraction) + (nParticles - 1.0);
}


void VariationalMonteCarlo::UpdateEnergies(const int &i)
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
        //waveFunctionNew = waveFunctionOld;  // Probably unnecessary since the wavefunction didn't change.
    }

    if (derivationType == "Analytical") {
        // Update energies (without numerical derivation and interaction)
        deltaEnergy = Hamiltonian::LocalEnergy(rNew, nParticles, nDimensions, alpha);
    } else if (derivationType == "Numerical")
    {
        // Update energies using numerical derivation
        deltaEnergy = Hamiltonian::NumericalLocalEnergy(rNew, nParticles, nDimensions, alpha, stepLength, beta);
    } else if (derivationType == "Interaction")
    {
        // Update energies using interaction
        deltaEnergy = Hamiltonian::LocalEnergyInteraction(rNew, nParticles, nDimensions, alpha, beta, a);
    }

    // Delete energies that should not have been accepted. Only affects interactions.
    if (fabs(deltaEnergy) > 100)
    {
        deltaEnergy        = 0;
        thrownCounter     += 1;
        acceptanceCounter -= 1;  // To account for thrown energies
    }

    energySum         += deltaEnergy;
    energySquaredSum  += (deltaEnergy*deltaEnergy);

    if (cycleType == "SteepestDescent")
    {
        deltaPsi           = Wavefunction::DerivativePsi(rNew, nParticles, nDimensions, beta);
        psiSum            += deltaPsi;
        psiTimesEnergySum += deltaPsi*deltaEnergy;
    }
}


double VariationalMonteCarlo::SteepestDescent(const int &nParticles)
{
    double eta                   = 0.001;
    int nAlpha                   = 500;
    double averagePsi            = 0.0;
    double averageEnergy         = 0.0;
    double averagePsiTimesEnergy = 0.0;
    double localEnergyDerivative = 0.0;
    double scalingFactor         = 1.0/(nCycles*nParticles);

    // Loop over number of alphas
    for (int i = 0; i < nAlpha; i++)
    {
        acceptanceCounter = 0;
        waveFunctionOld   = 0.0;
        waveFunctionNew   = 0.0;
        energySum         = 0.0;
        energySquaredSum  = 0.0;
        deltaEnergy       = 0.0;
        acceptanceWeight  = 0.0;
        psiSum            = 0.0;
        psiTimesEnergySum = 0.0;
        deltaPsi          = 0.0;

        // Run Monte Carlo cycles
        MonteCarloCycles();

        // Update energies:
        averagePsi    = psiSum*scalingFactor;
        averageEnergy = energySum*scalingFactor;
        averagePsiTimesEnergy = psiTimesEnergySum*scalingFactor;
        localEnergyDerivative = 2.0*(averagePsiTimesEnergy - averagePsi*averageEnergy);

        // Calculate alpha
        std::cout << "New alpha: "  << alpha
                  << " Average energy: " << averageEnergy << std::endl;
        std::cout << acceptanceCounter*scalingFactor << "   " << acceptanceCounter << std::endl;

        alpha -= eta * localEnergyDerivative;
    }
    return alpha;
}
