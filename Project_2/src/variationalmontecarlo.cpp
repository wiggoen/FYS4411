#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include "inc/hermite.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>    /* high resolution timing: http://en.cppreference.com/w/cpp/chrono/c/clock */
#include <string>
#include <cmath>
#include <stdlib.h>  /* Exit failure, to force the program to stop: exit(EXIT_FAILURE); */

#ifdef MPI_ON
#include <mpi.h>
#endif


VariationalMonteCarlo::VariationalMonteCarlo() :
    nDimensions(2)
{
    isOneBodySetup = false;
}


VariationalMonteCarlo::~VariationalMonteCarlo( void )
{

}


arma::rowvec VariationalMonteCarlo::RunVMC(const int nParticles, const int nCycles, const double alpha, const double beta,
                                           const double omega, const double stepLength, const double timeStep,
                                           const bool UseJastrowFactor, const bool UseImportanceSampling,
                                           const bool UseFermionInteraction, const bool UseAnalyticalExpressions,
                                           bool UseNumericalPotentialEnergy, std::string cycleType,
                                           const int cycleStepToFile, bool UseMPI, int terminalizationFactor)
{
    /* Adding variables to member variables */
    this->nParticles                  = nParticles;
    this->nCycles                     = nCycles;
    this->alpha                       = alpha;
    this->beta                        = beta;
    this->omega                       = omega;
    this->stepLength                  = stepLength;
    this->timeStep                    = timeStep;
    this->UseJastrowFactor            = UseJastrowFactor;
    this->UseImportanceSampling       = UseImportanceSampling;
    this->UseFermionInteraction       = UseFermionInteraction;
    this->UseAnalyticalExpressions    = UseAnalyticalExpressions;
    this->UseNumericalPotentialEnergy = UseNumericalPotentialEnergy;
    this->cycleType                   = cycleType;
    this->cycleStepToFile             = cycleStepToFile;
    this->terminalizationFactor       = terminalizationFactor;


    /* Initialize matrices and variables */
    rOld = arma::zeros<arma::mat>(nParticles, nDimensions);
    rNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    if (UseImportanceSampling)
    {
        QForceOld = arma::zeros<arma::mat>(nParticles, nDimensions);
        QForceNew = arma::zeros<arma::mat>(nParticles, nDimensions);
    }
    if (!UseAnalyticalExpressions)
    {
        waveFunctionOld = 0.0;
        waveFunctionNew = 0.0;
    }
    energySum             = 0.0;
    energySquaredSum      = 0.0;
    deltaEnergy           = 0.0;
    acceptanceWeight      = 0.0;
    acceptanceCounter     = 0.0;
    slaterRatio           = 0.0;
    jastrowRatio          = 0.0;
    energyVector          = arma::zeros<arma::rowvec>(nDimensions+1);
    numericalEnergyVector = arma::zeros<arma::rowvec>(nDimensions+1);
    kineticEnergySum      = 0.0;
    potentialEnergySum    = 0.0;


    /* Initial trial positions */
    if (!UseImportanceSampling) { InitialTrialPositionsBruteForce(rOld); }
    else                        { InitialTrialPositionsImportanceSampling(rOld); }
    rNew = rOld;

    /* Initial slater determinants */
    SlaterUpOld          = arma::zeros<arma::mat>(nParticles/2, nParticles/2);
    SlaterDownOld        = arma::zeros<arma::mat>(nParticles/2, nParticles/2);
    SlaterInitialization();
    SlaterUpNew          = SlaterUpOld;
    SlaterDownNew        = SlaterDownOld;
    InverseSlaterUpOld   = inv(SlaterUpOld);
    InverseSlaterDownOld = inv(SlaterDownOld);
    InverseSlaterUpNew   = InverseSlaterUpOld;
    InverseSlaterDownNew = InverseSlaterDownOld;

    if (UseImportanceSampling)
    {
        if (!UseAnalyticalExpressions)
        {
            /* Store the current value of the wave function and quantum force */
            waveFunctionOld = Wavefunction::NumericalTrialWaveFunction(rOld, nParticles, alpha, beta, omega, UseJastrowFactor);

            Wavefunction::NumericalQuantumForce(rOld, QForceOld, nParticles, nDimensions, alpha, beta, omega,
                                                UseJastrowFactor);
        } else
        {
            for (int i = 0; i < nParticles; i++)
            {
                /* Initialize quantum force */
                Wavefunction::QuantumForce(rOld, QForceOld, nParticles, nDimensions, alpha, beta, omega, UseJastrowFactor,
                                           InverseSlaterUpOld, InverseSlaterDownOld, i);
            }
        }
        QForceNew = QForceOld;
    }


    /* MPI */
#ifdef MPI_ON
    if (UseMPI)
    {
        /* Initialize the MPI environment */
        MPI_Init(NULL, NULL);

        /* Get number of processes */
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        /* Get the rank of the process */
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        /* Get the name of the processor */
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);

        /* Print off a hello world message */
        std::cout << "Hello world from processor " << processor_name << ", rank " << world_rank
                  << " out of " << world_size << " processors." << std::endl;
    }
#endif

    /* Start timing */
    auto start_time = std::chrono::high_resolution_clock::now();

    /* Run Monte Carlo cycles */
    if (cycleType == "MonteCarlo" || cycleType == "OneBodyDensity")
    {
        std::cout << "Running MC Cycles.." << std::endl;
        MonteCarloCycles();
    } else if (cycleType == "SteepestDescent") {
        SteepestDescent();
    }

    /* Timing finished */
    auto end_time = std::chrono::high_resolution_clock::now();
    double runTime = (std::chrono::duration<double> (end_time - start_time).count());

    /* Normalizing */
    double normalizationFactor = 1.0/(nCycles - terminalizationFactor);

    /* Calculation of averages */
    double energy = energySum * normalizationFactor;
    double energySquared = energySquaredSum * normalizationFactor;

    kineticEnergy   = kineticEnergySum * normalizationFactor;
    potentialEnergy = potentialEnergySum * normalizationFactor;

    double variance = (energySquared - energy*energy) * normalizationFactor;
    double acceptanceRatio = acceptanceCounter * normalizationFactor;


    if (cycleType == "OneBodyDensity")                                      /* SHOULD THERE BE A TERMINALIZATION FACTOR HERE??? */
    {
        for (int i = 0; i < nBins; i++)
        {
            hist(i) /= volume(i) * normalizationFactor;
        }
        std::ofstream histogram;  /* Write position matrix to file */
        histogram.open("../Project_2/results/histogram.txt");

        /* Loop over bins */
        for (int i = 0; i < nBins; i++)
        {
            /* Write to file */
            histogram << hist(i) << std::endl;
        }
        histogram.close();
    }

    /* Vector containing the results of the run */
    arma::rowvec runDetails;
    runDetails << runTime << energy << energySquared << variance << acceptanceRatio << kineticEnergy << potentialEnergy;

    /* Finalize the MPI environment */
#ifdef MPI_ON
    if (UseMPI) { MPI_Finalize(); }
#endif

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


void VariationalMonteCarlo::SlaterInitialization( void )
{
    const arma::mat QuantumNumber = Hermite::QuantumNumbers();

    /* initialize slater determinants */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            SlaterUpOld(i, j)   = Wavefunction::phi(rOld, alpha, omega, nx, ny, i);
            SlaterDownOld(i, j) = Wavefunction::phi(rOld, alpha, omega, nx, ny, i+nParticles/2);
        }
    }
}


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
    if (cycleStepToFile != 0 && world_rank == 0)
    {
        outputEnergy.open("../Project_2/results/energies.txt");
        outputDistance.open("../Project_2/results/distances.txt");
    }

    /* Loop over Monte Carlo cycles */
    for (int cycle = 0; cycle < nCycles; cycle++)
    {
        cycleNumber = cycle;

        /* Sampling */
        if (!UseImportanceSampling)
        {
            /* Brute force sampling */
            MetropolisBruteForce(rNew, rOld);
        } else
        {
            /* Importance sampling */
            ImportanceSampling(rNew, rOld, QForceNew, QForceOld);
        }

        /* Write to file */
        if (cycleStepToFile != 0 && world_rank == 0 && cycle % cycleStepToFile == 0)
        {
            if (cycle >= terminalizationFactor)
            {
                outputEnergy << std::setw(10) << cycle << "     " << std::setprecision(6) << deltaEnergy
                             << std::endl;

                if (nParticles == 2)
                {
                    double particleDistance = Hamiltonian::ParticleDistance(rNew.row(0), rNew.row(1));    /* need to fix this for more particles */
                    outputDistance << std::setw(10) << cycle << "     " << std::setprecision(6)
                                   << particleDistance << std::endl;
                }
            }
        }
    }
    /* Close output file */
    if (cycleStepToFile != 0 && world_rank == 0)
    {
        outputEnergy.close();
        outputDistance.close();
    }
}


int VariationalMonteCarlo::PickRandomParticle( void )
{
    static std::random_device rd;  /* Initialize the seed for the random number engine */
    static std::mt19937_64 gen(rd());  /* Call the Mersenne Twister algorithm */
    static std::uniform_int_distribution<int> dist_int(0, nParticles-1);
    return dist_int(gen);
}


void VariationalMonteCarlo::MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld)
{
    /* New position to test */
    int i = PickRandomParticle();
    for (int j = 0; j < nDimensions; j++)
    {
        rNew(i, j) = rOld(i, j) + (UniformRandomNumber() - 0.5) * stepLength;
    }

    /* Recalculate the value of the wave function ratio */
    if (!UseAnalyticalExpressions)
    {
        /* using numerical expressions */
        waveFunctionNew = Wavefunction::NumericalTrialWaveFunction(rNew, nParticles, alpha, beta, omega, UseJastrowFactor);

        acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
    } else
    {
        /* using analytical expressions */
        slaterRatio = Wavefunction::SlaterRatio(rNew, nParticles, alpha, omega, InverseSlaterUpOld, InverseSlaterDownOld, i);

        jastrowRatio = 1.0;
        if (UseJastrowFactor)
        {
            jastrowRatio = Wavefunction::JastrowRatio(rNew, rOld, nParticles, beta);
        }

        acceptanceWeight = slaterRatio*slaterRatio * jastrowRatio*jastrowRatio;
    }

    UpdateEnergies(i);
}


void VariationalMonteCarlo::ImportanceSampling(arma::mat &rNew, const arma::mat &rOld, arma::mat &QForceNew,
                                               const arma::mat &QForceOld)
{
    double diffusionCoefficient = 0.5;

    /* New position to test */
    int i = PickRandomParticle();
    for (int j = 0; j < nDimensions; j++)
    {
        rNew(i, j) = rOld(i, j) + diffusionCoefficient*QForceOld(i, j)*timeStep + GaussianRandomNumber()*sqrt(timeStep);
    }
    /* Recalculate the value of the wave function ratio and the quantum force */
    if (!UseAnalyticalExpressions)
    {
        /* using numerical expressions */
        waveFunctionNew = Wavefunction::NumericalTrialWaveFunction(rNew, nParticles, alpha, beta, omega, UseJastrowFactor);
        Wavefunction::NumericalQuantumForce(rNew, QForceNew, nParticles, nDimensions, alpha, beta, omega,
                                            UseJastrowFactor);

        double wavefunctionRatio = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
        double greensRatio = GreensRatio(rNew, rOld, QForceNew, QForceOld, diffusionCoefficient, timeStep, i);
        acceptanceWeight = greensRatio * wavefunctionRatio;
    } else
    {
        /* using analytical expressions */
        Wavefunction::QuantumForce(rNew, QForceNew, nParticles, nDimensions, alpha, beta, omega, UseJastrowFactor,
                                   InverseSlaterUpOld, InverseSlaterDownOld, i);

        slaterRatio = Wavefunction::SlaterRatio(rNew, nParticles, alpha, omega, InverseSlaterUpOld, InverseSlaterDownOld, i);

        jastrowRatio = 1.0;
        if (UseJastrowFactor)
        {
            jastrowRatio = Wavefunction::JastrowRatio(rNew, rOld, nParticles, beta);
        }
        double greensRatio = GreensRatio(rNew, rOld, QForceNew, QForceOld, diffusionCoefficient, timeStep, i);
        acceptanceWeight = slaterRatio*slaterRatio * jastrowRatio*jastrowRatio * greensRatio;
    }

    UpdateEnergies(i);
}


double VariationalMonteCarlo::GreensRatio(const arma::mat &rNew, const arma::mat &rOld, const arma::mat &QForceNew,
                                          const arma::mat &QForceOld, const double &diffusionCoefficient,
                                          const double &timeStep, const int &i)
{
    double fraction  = 1.0/(4.0*diffusionCoefficient*timeStep);
    arma::rowvec yx  = rNew.row(i) - rOld.row(i) - diffusionCoefficient*timeStep*QForceOld.row(i);
    double yxSquared = arma::dot(yx, yx);
    arma::rowvec xy  = rOld.row(i) - rNew.row(i) - diffusionCoefficient*timeStep*QForceNew.row(i);
    double xySquared = arma::dot(xy, xy);
    return exp((-xySquared + yxSquared) * fraction) + (nParticles - 1.0);
}


void VariationalMonteCarlo::UpdateInverseSlater(const int &i)
{
    const arma::mat QuantumNumber = Hermite::QuantumNumbers();

    double factor = 1.0/slaterRatio;

    if (i < nParticles/2)
    {
        /* Update InverseSlaterUp */
        for (int k = 0; k < nParticles/2; k++)
        {
            /* i-th column */
            InverseSlaterUpNew(k, i) = factor*InverseSlaterUpOld(k, i);
        }

        for (int j = 0; j < nParticles/2; j++)
        {
            /* j-th column */
            if (j != i)
            {
                double Sj = 0.0;

                for (int l = 0; l < nParticles/2; l++)
                {
                    int nx = QuantumNumber(l, 0);
                    int ny = QuantumNumber(l, 1);
                    Sj += Wavefunction::phi(rNew, alpha, omega, nx, ny, i)*InverseSlaterUpOld(l, j);
                }

                for (int k = 0; k < nParticles/2; k++)
                {
                    InverseSlaterUpNew(k, j) = InverseSlaterUpOld(k, j) - (Sj/slaterRatio)*InverseSlaterUpOld(k, i);
                }
            }
        }
    } else
    {
        /* Update InverseSlaterDown */
        for (int k = 0; k < nParticles/2; k++)
        {
            /* i-th column */
            InverseSlaterDownNew(k, i-nParticles/2) = factor*InverseSlaterDownOld(k, i-nParticles/2);
        }

        for (int j = 0; j < nParticles/2; j++)
        {
            /* j-th column */
            if (j != (i-nParticles/2))
            {
                double Sj = 0.0;

                for (int l = 0; l < nParticles/2; l++)
                {
                    int nx = QuantumNumber(l, 0);
                    int ny = QuantumNumber(l, 1);
                    Sj += Wavefunction::phi(rNew, alpha, omega, nx, ny, i)*InverseSlaterDownOld(l, j);
                }

                for (int k = 0; k < nParticles/2; k++)
                {
                    InverseSlaterDownNew(k, j) = InverseSlaterDownOld(k, j) - (Sj/slaterRatio)*InverseSlaterDownOld(k, i-nParticles/2);
                }
            }
        }
    }
    InverseSlaterUpOld   = InverseSlaterUpNew;
    InverseSlaterDownOld = InverseSlaterDownNew;
}


void VariationalMonteCarlo::UpdateEnergies(const int &i)
{
    /* Test is performed by moving one particle at the time. Accept or reject this move. */
    if (UniformRandomNumber() <= acceptanceWeight)
    {
        rOld.row(i) = rNew.row(i);
        if (UseImportanceSampling)    { QForceOld.row(i) = QForceNew.row(i); }
        if (UseAnalyticalExpressions) { UpdateInverseSlater(i); }
        else                          { waveFunctionOld = waveFunctionNew; }
        if (cycleNumber >= terminalizationFactor) { acceptanceCounter += 1; }
    } else
    {
        rNew.row(i) = rOld.row(i);
        if (UseImportanceSampling)    { QForceNew.row(i) = QForceOld.row(i); }
        if (UseAnalyticalExpressions)
        {
            if (i < nParticles/2) { InverseSlaterUpNew.col(i)                = InverseSlaterUpOld.col(i); }
            else                  { InverseSlaterDownNew.col(i-nParticles/2) = InverseSlaterDownOld.col(i-nParticles/2); }
        }
    }

    /* Update energies */
    if (!UseAnalyticalExpressions)
    {
        /* using numerical expressions */
        numericalEnergyVector = Hamiltonian::NumericalLocalEnergy(rNew, nParticles, nDimensions, alpha, beta, omega,
                                                                  UseJastrowFactor, UseFermionInteraction,
                                                                  UseNumericalPotentialEnergy);
        if (cycleNumber >= terminalizationFactor)
        {
            deltaEnergy         = numericalEnergyVector.at(0);
            kineticEnergySum   += numericalEnergyVector.at(1);
            potentialEnergySum += numericalEnergyVector.at(2);
        }
    } else
    {
        /* using analytical expressions */
        energyVector = Hamiltonian::LocalEnergy(rNew, nParticles, nDimensions, alpha, beta, omega,
                                                UseJastrowFactor, UseFermionInteraction,
                                                InverseSlaterUpNew, InverseSlaterDownNew);
        if (cycleNumber >= terminalizationFactor)
        {
            deltaEnergy         = energyVector.at(0);
            kineticEnergySum   += energyVector.at(1);
            potentialEnergySum += energyVector.at(2);
        }
    }
    energySum        += deltaEnergy;
    energySquaredSum += (deltaEnergy*deltaEnergy);

    if (cycleType == "SteepestDescent" && cycleNumber >= terminalizationFactor)
    {
        if (nParticles == 2)
        {
            dPsiOfAlpha  = Wavefunction::DerivativePsiOfAlpha(rNew, omega);
            dPsiOfBeta   = Wavefunction::DerivativePsiOfBeta(rNew, beta);
        } else
        {
            dPsiOfAlpha  = Wavefunction::DerivativePsiManyOfAlpha(rNew, nParticles, alpha, omega);
            dPsiOfBeta   = Wavefunction::DerivativePsiManyOfBeta(rNew, nParticles, beta);
        }

        psiSumAlpha += dPsiOfAlpha;
        psiSumBeta  += dPsiOfBeta;
        psiOfAlphaTimesEnergySum += dPsiOfAlpha*deltaEnergy;
        psiOfBetaTimesEnergySum  += dPsiOfBeta*deltaEnergy;
    }

    if (cycleType == "OneBodyDensity")
    {
        OneBodyDensity();
    }
}


void VariationalMonteCarlo::SteepestDescent( void )
{
    std::cout << "Running steepest descent..." << std::endl;
    double eta                   = 0.001;
    int nAlpha                   = 1000;
    double averagePsiAlpha       = 0.0;
    double averagePsiBeta        = 0.0;
    double averageEnergy         = 0.0;
    double alphaEnergyDerivative = 0.0;
    double betaEnergyDerivative  = 0.0;
    double scalingFactor         = 1.0/(nCycles - terminalizationFactor);
    double averagePsiOfAlphaTimesEnergy = 0.0;
    double averagePsiOfBetaTimesEnergy  = 0.0;

    std::cout << "Iteration " << " New_alpha " << " New_beta " << "   Energy " << std::endl;
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
        psiOfAlphaTimesEnergySum = 0.0;
        psiOfBetaTimesEnergySum  = 0.0;

        /* Run Monte Carlo cycles */
        MonteCarloCycles();

        /* Update energies */
        averagePsiAlpha = psiSumAlpha*scalingFactor;
        averagePsiBeta  = psiSumBeta*scalingFactor;
        averageEnergy   = energySum*scalingFactor;
        averagePsiOfAlphaTimesEnergy = psiOfAlphaTimesEnergySum*scalingFactor;
        averagePsiOfBetaTimesEnergy  = psiOfBetaTimesEnergySum*scalingFactor;
        alphaEnergyDerivative = 2.0*(averagePsiOfAlphaTimesEnergy - averagePsiAlpha*averageEnergy);
        betaEnergyDerivative  = 2.0*(averagePsiOfBetaTimesEnergy - averagePsiBeta*averageEnergy);

        if (alphaEnergyDerivative == 0 && betaEnergyDerivative == 0 && i>10)
        {
            std::cout << "Optimal parameters found, stopping steepest descent" << std::endl;
            i=nAlpha-1;
        }

        /* Calculate alpha and beta */
        alpha -= eta * alphaEnergyDerivative;
        beta  -= eta * betaEnergyDerivative;

        std::cout << std::setw(9) << i << std::setw(11) << alpha << std::setw(10) << beta
                  << std::setw(10) << averageEnergy << std::endl;
    }
}


void VariationalMonteCarlo::SetOneBody( void )
{
    int max_r = 10;
    nBins = 800;
    r_step = (double) max_r/nBins;
    hist = arma::zeros<arma::rowvec>(nBins);
    volume = arma::zeros<arma::rowvec>(nBins);

    volume(0) = r_step*r_step;
    for (int i = 1; i < nBins; i++)
    {
        volume(i) = pow(r_step*(i+1), 2);
    }
}


void VariationalMonteCarlo::OneBodyDensity( void )
{
    if (!isOneBodySetup)
    {
        SetOneBody();
        isOneBodySetup = true;
    }

    for (int i = 0; i < nParticles; i++)
    {
        double rNorm = arma::norm(rNew.row(i));
        for (int j = 1; j < nBins; j++)
        {
            if ( (rNorm < j*r_step) && (rNorm >= (r_step*(j-1))) )
            {
                hist(j-1) += 1;
            }
        }
    }

}
