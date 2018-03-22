#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>  // high resolution timing: http://en.cppreference.com/w/cpp/chrono/c/clock
#include <cmath>
//#include <stdlib.h> /* Exit failure <- to force the program to stop: exit(EXIT_FAILURE);  */


VariationalMonteCarlo::VariationalMonteCarlo()
{

}


VariationalMonteCarlo::~VariationalMonteCarlo()
{

}


arma::rowvec VariationalMonteCarlo::RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles,
                                                             double alpha, double stepLength, double timeStep,
                                                             int cycleStepToFile)
{
    // Adding variables to member variables
    this->nParticles = nParticles;
    this->nDimensions = nDimensions;
    this->nCycles = nCycles;
    this->alpha = alpha;
    this->stepLength = stepLength;
    this->timeStep = timeStep;
    this->cycleStepToFile = cycleStepToFile;

    // For writing to file										<<<< REMEMBER TO REMOVE THESE!
	/*
	if ( nParticles == 1 ) { strParticles = "1"; }
    if ( nParticles == 10 ) { strParticles = "10"; }
    if ( nParticles == 100 ) { strParticles = "100"; }
    if ( nParticles == 500 ) { strParticles = "500"; }
	if ( nDimensions == 1 ) { strDimensions = "1"; }
	if ( nDimensions == 2 ) { strDimensions = "2"; }
    if ( nDimensions == 3 ) { strDimensions = "3"; }
    if ( alpha == 0.1 ) { strAlpha = "01"; }
    if ( alpha == 0.2 ) { strAlpha = "02"; }
    if ( alpha == 0.3 ) { strAlpha = "03"; }
    if ( alpha == 0.4 ) { strAlpha = "04"; }
    if ( alpha == 0.5 ) { strAlpha = "05"; }
    if ( alpha == 0.6 ) { strAlpha = "06"; }
    if ( alpha == 0.7 ) { strAlpha = "07"; }
    if ( alpha == 0.8 ) { strAlpha = "08"; }
    if ( alpha == 0.9 ) { strAlpha = "09"; }
	*/

    // Initialize matrices and variables
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    waveFunctionOld = 0;
    waveFunctionNew = 0;
    energySum = 0;
    energySquaredSum = 0;
    deltaEnergy = 0;
    acceptanceWeight = 0;
    acceptanceCounter = 0;
    psiSum = 0;
    psiTimesEnergySum = 0;
    //psiSum = 0;

    arma::rowvec runDetails;
    double runTime = 0;
    double energy = 0;
    double energySquared = 0;
    double variance = 0;
    double acceptanceRatio = 0;
    double normalizationFactor = 1.0/(nCycles * nParticles);

    // CHOOSE SAMPLING METHOD                    <<< ---
    samplingType = "BruteForce";
    //samplingType = "Importance";

    // Initial trial positions
    if (samplingType == "BruteForce") {
        InitialTrialPositionsBruteForce(rOld);
    } else if (samplingType == "Importance")
    {
        InitialTrialPositionsImportanceSampling(rOld);
    }
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
    runTime = (std::chrono::duration<double> (end_time - start_time).count());

    // Calculation
    energy = energySum * normalizationFactor;
    energySquared = energySquaredSum * normalizationFactor;

    variance = (energySquared - energy*energy) * normalizationFactor;
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


void VariationalMonteCarlo::MonteCarloCycles()
{
    std::ofstream myfile;
    myfile.open("../Project_1/results.txt");
    //myfile.open("../Project_1/results/results_" + strParticles + "p_" + strDimensions + "d_alpha" + strAlpha + ".txt");

    // Loop over Monte Carlo cycles
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        // Sampling
        if (samplingType == "BruteForce") {
            MetropolisBruteForce(rNew, rOld, waveFunctionOld, waveFunctionNew);
        } else if (samplingType == "Importance")
        {
            ImportanceSampling(rNew, rOld, QForceOld, QForceNew, waveFunctionOld, waveFunctionNew);
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
        // Analytic:
        //waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, nParticles, nDimensions, alpha);
        // Numeric:


        acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);

        UpdateEnergies(i);
    }
}

double VariationalMonteCarlo::NumericalDerivation(arma::mat R, double &waveFunctionNew)
{
    waveFunctionNew = Wavefunction::TrialWaveFunction(R, nParticles, nDimensions, alpha);
    double Rsquared = 0;
    double kineticEnergy = 0;
    double externalPotential =0;
    double h = stepLength;
    double h2 = h*h;

    // Kinetic energy
    arma::mat Rplus = R + h;
    arma::mat Rminus = R - h;

    for (int i = 0; i < nParticles; i++)
    {
        Rsquared = 0;
        for (int j = 0; j < nParticles; j++)
        {
            Rsquared += R(i,j)*R(i,j);
        }
        externalPotential += 0.5*Rsquared;
    }
    double waveFunctionPlus = Wavefunction::TrialWaveFunction(Rplus, nParticles, nDimensions, alpha);
    double waveFunctionMinus = Wavefunction::TrialWaveFunction(Rminus, nParticles, nDimensions, alpha);
    kineticEnergy -= (waveFunctionPlus+waveFunctionMinus - 2*waveFunctionNew);
    kineticEnergy = 0.5*kineticEnergy*h2/waveFunctionNew;
    return kineticEnergy + externalPotential;
}



void VariationalMonteCarlo::ImportanceSampling(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld,
                                               arma::mat &QForceNew, double &waveFunctionOld,
                                               double &waveFunctionNew)
{
    double diffusionCoefficient = 0.5;
    double wavefunctionsSquared = 0;
    double GreensRatio = 0;

    // New position to test
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rNew(i, j) = rOld(i, j) + diffusionCoefficient*QForceOld(i, j)*timeStep + GaussianRandomNumber()*sqrt(timeStep);
        }

        // Recalculate the value of the wave function and the quantum force
        waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, nParticles, nDimensions, alpha);
        Wavefunction::QuantumForce(rNew, QForceNew, alpha);

        wavefunctionsSquared = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
        GreensRatio = GreensFunction(rOld, rNew, QForceOld, diffusionCoefficient, timeStep, i);
        acceptanceWeight = GreensRatio * wavefunctionsSquared;

        UpdateEnergies(i);
    }
}




double VariationalMonteCarlo::GreensFunction(const arma::mat &rOld, const arma::mat &rNew, const arma::mat &QForceOld,
                                             double &diffusionCoefficient, double &timeStep, int &i)
{
    double fraction = 1.0/(4.0*diffusionCoefficient*timeStep);
    arma::rowvec yx = rNew.row(i) - rOld.row(i) - diffusionCoefficient*timeStep*QForceOld.row(i);
    double yxSquared = arma::dot(yx, yx);
    arma::rowvec xy = rOld.row(i) - rNew.row(i) - diffusionCoefficient*timeStep*QForceNew.row(i);
    double xySquared = arma::dot(xy, xy);
    return exp((-xySquared + yxSquared) * fraction) + (nParticles - 1.0);
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
        waveFunctionNew = waveFunctionOld;              // SIKKERT UNØDVENDIG
    }

    // Update energies
    deltaEnergy        = Hamiltonian::LocalEnergy(rNew, nParticles, nDimensions, alpha);
    energySum         += deltaEnergy;
    energySquaredSum  += deltaEnergy*deltaEnergy;
    psiSum            += waveFunctionNew;
    psiTimesEnergySum += waveFunctionNew*deltaEnergy;
    //psiEnergyTot         = waveFunctionNew*energySum*scalingFactor;
    //averageEnergy        = totalEnergy;
    //averageEnergySquared = totalEnergySquared;

    //double psiEnergyAverage = psiEnergyTot;
    //double averagePsi = totalPsi;
    //double localEnergyDerivative = 2*(psiEnergyAverage - averagePsi*averageEnergy);
}


double VariationalMonteCarlo::SteepestDescent(int nParticles, int nDimensions, double initialAlpha)
{
    double eta              = 0.001;
    double nAlpha           = 50;
    double totalEnergy      = 0;
    double totalPsi         = 0;
    double psiEnergyTot     = 0;
    double averageEnergy    = 0;
    double averageEnergySquared;
    //arma::vec R = rOld;

    double scalingFactor = 1.0/(nParticles*nDimensions*nCycles);

    // Loop over number of alphas
    for (int i = 0; i<nAlpha; i++)
    {
        std::cout << energySum*scalingFactor << std::endl;
        // Run Monte Carlo cycles
        RunMonteCarloIntegration(nParticles, nDimensions, nCycles, alpha, stepLength, timeStep, cycleStepToFile);

        std::cout << energySum*scalingFactor << std::endl;

        // Update energies:
        totalEnergy          = energySum*scalingFactor;
        //totalEnergySquared   = energySquaredSum*scalingFactor;
        totalPsi             = psiSum*scalingFactor;
        psiEnergyTot         = psiTimesEnergySum*scalingFactor;
        averageEnergy        = totalEnergy;

        double psiEnergyAverage = psiEnergyTot;
        double averagePsi = totalPsi;
        double localEnergyDerivative = 2*(psiEnergyAverage - averagePsi*averageEnergy);

        // Print energies to check values:
        std::cout << "Average psi*energy: "      << psiEnergyAverage      << std::endl
                  << "Average psi:        "      << averagePsi            << std::endl
                  << "Average energy:     "      << averageEnergy         << std::endl
                  << "Local energy derivative: " << localEnergyDerivative << std::endl;

        alpha = alpha - eta*localEnergyDerivative;
        std::cout << alpha - eta*localEnergyDerivative << std::endl;
        std::cout << "New alpha: " << alpha << std::endl;
    }

    std::cout << std::fixed << std::setprecision(15) << alpha << std::endl;
    return alpha;
}
