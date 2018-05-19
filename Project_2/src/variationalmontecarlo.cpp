#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>    /* high resolution timing: http://en.cppreference.com/w/cpp/chrono/c/clock */
#include <string>
#include <cmath>
#include <stdlib.h>  /* Exit failure, to force the program to stop: exit(EXIT_FAILURE); */
#include <mpi.h>
//#include <stdio.h>


VariationalMonteCarlo::VariationalMonteCarlo() :
    nDimensions(2)
{

}


VariationalMonteCarlo::~VariationalMonteCarlo( void )
{

}


arma::rowvec VariationalMonteCarlo::RunVMC(const int nParticles, const int nCycles, const double alpha, const double beta,
                                           const double omega, const double spinParameter, const double stepLength,
                                           const double timeStep, const bool UseJastrowFactor,
                                           const bool UseImportanceSampling, const bool UseFermionInteraction,
                                           const bool UseAnalyticalExpressions, bool UseNumericalPotentialEnergy,
                                           std::string cycleType, const int cycleStepToFile, bool UseMPI)
{
    /* Adding variables to member variables */
    this->nParticles                  = nParticles;
    this->nCycles                     = nCycles;
    this->alpha                       = alpha;
    this->beta                        = beta;
    this->omega                       = omega;
    this->spinParameter               = spinParameter;
    this->stepLength                  = stepLength;
    this->timeStep                    = timeStep;
    this->UseJastrowFactor            = UseJastrowFactor;
    this->UseImportanceSampling       = UseImportanceSampling;
    this->UseFermionInteraction       = UseFermionInteraction;
    this->UseAnalyticalExpressions    = UseAnalyticalExpressions;
    this->UseNumericalPotentialEnergy = UseNumericalPotentialEnergy;
    this->cycleType                   = cycleType;
    this->cycleStepToFile             = cycleStepToFile;


    /* Initialize matrices and variables */
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

    arma::rowvec runDetails;
    double runTime = 0.0;
    double energy = 0.0;
    double energySquared = 0.0;
    double variance = 0.0;
    double acceptanceRatio = 0.0;


    /* Initial trial positions */
    if (!UseImportanceSampling) { InitialTrialPositionsBruteForce(rOld); }
    else                        { InitialTrialPositionsImportanceSampling(rOld); }
    rNew = rOld;

    //if (UseFermionInteraction)  { CheckInitialDistance(rOld); }

    /* Store the current value of the wave function and quantum force */
    waveFunctionOld = Wavefunction::TrialWaveFunction(rOld, alpha, beta, omega, spinParameter, UseJastrowFactor);
    Wavefunction::QuantumForce(rOld, QForceOld, alpha, beta, omega, spinParameter, UseJastrowFactor);
    QForceNew = QForceOld;

    /* MPI */
    if (UseMPI)
    {
        /* Initialize the MPI environment */
        MPI_Init(NULL, NULL);

        /* Get number of processes */
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        /* Get the rank of the process */
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        /* Get the name of the processor */
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);

        /* Print off a hello world message */
        //printf("Hello world from processor %s, rank %d"
        //       " out of %d processors\n",
        //       processor_name, world_rank, world_size);
        std::cout << "Hello world from processor " << processor_name << ", rank " << world_rank
                  << " out of " << world_size << " processors." << std::endl;
    }

    /* Start timing */
    auto start_time = std::chrono::high_resolution_clock::now();

    /* Run Monte Carlo cycles */
    if (cycleType == "MonteCarlo")//|| cycleType == "OneBodyDensity")
    {
        std::cout << "Running MC Cycles.." << std::endl;
        MonteCarloCycles();
    } else if (cycleType == "SteepestDescent") {
        SteepestDescent(nParticles);
    }

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


    /* Vector containing the results of the run */
    runDetails << runTime << energy << energySquared << variance << acceptanceRatio;

    /* Finalize the MPI environment */
    if (UseMPI) { MPI_Finalize(); }

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


//void VariationalMonteCarlo::RedrawPosition(arma::mat &r, const int &i)
//{
//    for (int j = 0; j < nDimensions; j++)
//    {
//        if (!UseImportanceSampling)
//        {
//            /* Brute force sampling */
//            r(i, j) = (UniformRandomNumber() - 0.5) * stepLength;
//        } else
//        {
//            /* Importance sampling */
//            r(i, j) = GaussianRandomNumber()*sqrt(timeStep);
//        }
//    }
//}


//void VariationalMonteCarlo::CheckInitialDistance(arma::mat &rOld)
//{
//    double distance = 0;
//    for (int i = 0; i < nParticles; i++)
//    {
//        for (int j = i+1; j < nParticles; j++)
//        {
//            distance = Hamiltonian::ParticleDistance(rOld.row(i), rOld.row(j));
//            while (distance < a)
//            {
//                /* Making sure that no particle lies within a distance 'a' from each other */
//                RedrawPositionImportanceSampling(rOld, i);
//                distance = Hamiltonian::ParticleDistance(rOld.row(i), rOld.row(j));
//            }
//        }
//    }
//}


double VariationalMonteCarlo::UniformRandomNumber( void )
{
    static std::random_device rd;  /* Initialize the seed for the random number engine */
    static std::mt19937_64 gen(rd());  /* Call the Mersenne Twister algorithm */
    /* Set up the uniform distribution for x in [0, 1] */
    static std::uniform_real_distribution<double> UniformNumberGenerator(0.0,1.0);
    return UniformNumberGenerator(gen);
}


double VariationalMonteCarlo::GaussianRandomNumber( void )
{
    static std::random_device rd;  /* Initialize the seed for the random number engine */
    static std::mt19937_64 gen(rd());  /* Call the Mersenne Twister algorithm */
    /* Set up the normal distribution for x in [0, 1] */
    static std::normal_distribution<double> NormalDistribution(0.0,1.0);
    return NormalDistribution(gen);
}


void VariationalMonteCarlo::MonteCarloCycles( void )
{
    /* Open outputfiles to write */
    std::ofstream outputEnergy;
    std::ofstream outputDistance;
    if (cycleStepToFile != 0)
    {
        outputEnergy.open("../Project_2/energies.txt");
        outputDistance.open("../Project_2/distances.txt");
    }

    /* Loop over Monte Carlo cycles */
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        /* Sampling */
        if (!UseImportanceSampling)
        {
            /* Brute force sampling */
            MetropolisBruteForce(rNew, rOld, waveFunctionNew, waveFunctionOld);
        } else
        {
            /* Importance sampling */
            ImportanceSampling(rNew, rOld, QForceNew, QForceOld, waveFunctionNew, waveFunctionOld);
        }

        /* Write to file */
        if (cycleStepToFile != 0 && cycle % cycleStepToFile == 0)
        {
            if (cycle > 0)
            {
                outputEnergy << std::setw(10) << cycle << "     " << std::setprecision(6) << energySum/(cycle*nParticles)
                             << std::endl;
            }

            double particleDistance = Hamiltonian::ParticleDistance(rNew.row(0), rNew.row(1));
            outputDistance << std::setw(10) << cycle << "     " << std::setprecision(6)
                           << particleDistance << std::endl;
        }
    }
    /* Close output file */
    if (cycleStepToFile != 0)
    {
        outputEnergy.close();
        outputDistance.close();
    }
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
        waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, alpha, beta, omega, spinParameter, UseJastrowFactor);

        acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
        //std::cout << "wOld = " << waveFunctionOld << std::endl;
        //std::cout << "wNew = " <<waveFunctionNew << std::endl;
        //std::cout << "acceptW = " <<acceptanceWeight << std::endl;

        UpdateEnergies(i);
    }
}

void VariationalMonteCarlo::ImportanceSampling(arma::mat &rNew, const arma::mat &rOld, arma::mat &QForceNew,
                                               const arma::mat &QForceOld, double &waveFunctionNew,
                                               const double &waveFunctionOld)
{
    double diffusionCoefficient = 0.5;

    /* New position to test */
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rNew(i, j) = rOld(i, j) + diffusionCoefficient*QForceOld(i, j)*timeStep + GaussianRandomNumber()*sqrt(timeStep);
        }
        /* Recalculate the value of the wave function and the quantum force */
        waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, alpha, beta, omega, spinParameter, UseJastrowFactor);
        if (!UseAnalyticalExpressions)
        {
            /* using numerical expressions */
            Wavefunction::NumericalQuantumForce(rNew, QForceNew, nParticles, nDimensions, alpha, beta, omega, spinParameter,
                                                UseJastrowFactor);
        } else
        {
            /* using analytical expressions */
            Wavefunction::QuantumForce(rNew, QForceNew, alpha, beta, omega, spinParameter, UseJastrowFactor);
        }
        double wavefunctionRatio = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
        double greensRatio = GreensRatio(rNew, rOld, QForceNew, QForceOld, diffusionCoefficient, timeStep, i);
        acceptanceWeight = greensRatio * wavefunctionRatio;

        UpdateEnergies(i);
    }
}


double VariationalMonteCarlo::GreensRatio(const arma::mat &rNew, const arma::mat &rOld, const arma::mat &QForceNew,
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
    /* Test is performed by moving one particle at the time. Accept or reject this move. */
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
        //waveFunctionNew = waveFunctionOld;  /* Probably unnecessary since the wavefunction didn't change. */
    }

    /* Update energies */
    if (!UseAnalyticalExpressions)
    {
        /* using numerical expressions */
        deltaEnergy = Hamiltonian::NumericalLocalEnergy(rNew, nParticles, nDimensions, alpha, beta, omega, spinParameter,
                                                        UseJastrowFactor, UseNumericalPotentialEnergy);
    } else
    {
        /* using analytical expressions */
        deltaEnergy = Hamiltonian::LocalEnergy(rNew, nParticles, alpha, beta, omega, spinParameter, UseJastrowFactor,
                                               UseFermionInteraction);
    }
    energySum         += deltaEnergy;
    energySquaredSum  += (deltaEnergy*deltaEnergy);

    if (cycleType == "SteepestDescent")
    {
        dPsiOfAlpha        = Wavefunction::DerivativePsiOfAlpha(rNew, alpha, omega);
        dPsiOfBeta         = Wavefunction::DerivativePsiOfBeta(rNew, alpha, omega);
        psiSumAlpha       += dPsiOfAlpha;
        psiSumBeta        += dPsiOfBeta;
        psiOfAlphaTimesEnergySum   += dPsiOfAlpha*deltaEnergy;
        psiOfBetaTimesEnergySum   += dPsiOfBeta*deltaEnergy;
    }

    /*
    if (cycleType == "OneBodyDensity")
    {
        OneBodyDensity();
    }*/
}


double VariationalMonteCarlo::SteepestDescent(const int &nParticles)
{
    std::cout << "Running steepest descent..." << std::endl;
    double eta                   = 0.1;
    int nAlpha                   = 1000;
    double averagePsiAlpha       = 0.0;
    double averagePsiBeta        = 0.0;
    double averageEnergy         = 0.0;
    double alphaEnergyDerivative = 0.0;
    double betaEnergyDerivative  = 0.0;
    double scalingFactor         = 1.0/(nCycles*nParticles);
    double averagePsiOfAlphaTimesEnergy = 0.0;
    double averagePsiOfBetaTimesEnergy = 0.0;

    std::cout << "New alpha:    New beta:" << std::endl;
    /* Loop over number of alphas */
    for (int i = 0; i < nAlpha; i++)
    {
        acceptanceCounter = 0;
        waveFunctionOld   = 0.0;
        waveFunctionNew   = 0.0;
        energySum         = 0.0;
        energySquaredSum  = 0.0;
        deltaEnergy       = 0.0;
        acceptanceWeight  = 0.0;
        psiSumAlpha       = 0.0;
        psiSumBeta        = 0.0;
        deltaPsi          = 0.0;
        psiOfAlphaTimesEnergySum = 0.0;
        psiOfBetaTimesEnergySum  = 0.0;


        /* Run Monte Carlo cycles */
        MonteCarloCycles();

        /* Update energies */
        averagePsiAlpha    = psiSumAlpha*scalingFactor;
        averagePsiBeta     = psiSumBeta*scalingFactor;
        averageEnergy      = energySum*scalingFactor;
        averagePsiOfAlphaTimesEnergy = psiOfAlphaTimesEnergySum*scalingFactor;
        averagePsiOfBetaTimesEnergy  = psiOfBetaTimesEnergySum*scalingFactor;
        alphaEnergyDerivative = 2.0*(averagePsiOfAlphaTimesEnergy - averagePsiAlpha*averageEnergy);
        betaEnergyDerivative  = 2.0*(averagePsiOfBetaTimesEnergy - averagePsiBeta*averageEnergy);

        if (alphaEnergyDerivative == 0 && betaEnergyDerivative ==0)
        {
            std::cout << "Optimal parameters found, stopping steepest descent" << std::endl;
            i=nAlpha-1;
        }

        //std::cout << averagePsiOfAlphaTimesEnergy << std::endl;
        /* Calculate alpha and beta */
        alpha -= eta * alphaEnergyDerivative;
        beta  -= eta * betaEnergyDerivative;
        std::cout << alpha << "     " << beta  << std::endl;
    }
    return alpha;
}
