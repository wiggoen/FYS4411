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
    slaterRatio       = 0.0;
    jastrowRatio      = 0.0;
    if (nParticles > 2 || !UseAnalyticalExpressions)
    {
        energyVector          = arma::zeros<arma::rowvec>(nDimensions+1);
        numericalEnergyVector = arma::zeros<arma::rowvec>(nDimensions+1);
        kineticEnergySum   = 0.0;
        potentialEnergySum = 0.0;
    }

    /* Initial trial positions */
    if (!UseImportanceSampling) { InitialTrialPositionsBruteForce(rOld); }
    else                        { InitialTrialPositionsImportanceSampling(rOld); }
    rNew = rOld;

    /* Initial slater determinants and quantum numbers */
    if (nParticles > 2)
    {
        SlaterUpOld          = arma::zeros<arma::mat>(nParticles/2, nParticles/2);
        SlaterDownOld        = arma::zeros<arma::mat>(nParticles/2, nParticles/2);
        spinMatrix           = arma::zeros<arma::mat>(nParticles, nParticles);
        QuantumNumber        = Hermite::QuantumNumbers();
        SlaterInitialization();
        SlaterUpNew          = SlaterUpOld;
        SlaterDownNew        = SlaterDownOld;
        InverseSlaterUpOld   = inv(SlaterUpOld);
        InverseSlaterDownOld = inv(SlaterDownOld);
        InverseSlaterUpNew   = InverseSlaterUpOld;
        InverseSlaterDownNew = InverseSlaterDownOld;
    }

    /* Store the current value of the wave function and quantum force */
    if (UseImportanceSampling)
    {
        waveFunctionOld = Wavefunction::TrialWaveFunction(rOld, alpha, beta, omega, spinParameter, UseJastrowFactor);
        Wavefunction::QuantumForce(rOld, QForceOld, alpha, beta, omega, spinParameter, UseJastrowFactor);
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
        SteepestDescent(nParticles);
    }

    /* Timing finished */
    auto end_time = std::chrono::high_resolution_clock::now();
    double runTime = (std::chrono::duration<double> (end_time - start_time).count());

    /* Normalizing */
    double normalizationFactor = 1.0/(nCycles * nParticles);

    /* Calculation of averages */
    double energy = energySum * normalizationFactor;
    double energySquared = energySquaredSum * normalizationFactor;

    if (nParticles > 2 || !UseAnalyticalExpressions)
    {
        kineticEnergy   = kineticEnergySum * normalizationFactor;
        potentialEnergy = potentialEnergySum * normalizationFactor;
    }

    double variance = (energySquared - energy*energy);
    double acceptanceRatio = acceptanceCounter * normalizationFactor;


    if (cycleType == "OneBodyDensity")
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
    if (nParticles > 2 || !UseAnalyticalExpressions)
    {
        runDetails << runTime << energy << energySquared << variance << acceptanceRatio << kineticEnergy << potentialEnergy;
    } else
    {
        runDetails << runTime << energy << energySquared << variance << acceptanceRatio;
    }

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

    /* set up spin parameter matrix */
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            if (i < nParticles/2) { if (j < nParticles/2) { spinMatrix(i, j) = 1.0/3.0; }
                                    else                  { spinMatrix(i, j) = 1.0;     } }
            else                  { if (j < nParticles/2) { spinMatrix(i, j) = 1.0;     }
                                    else                  { spinMatrix(i, j) = 1.0/3.0; } }
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
        /* Sampling */
        if (nParticles > 2)
        {
            /* Only implemented brute-force Metropolis for many particles  */
            MetropolisBruteForce(rNew, rOld, waveFunctionNew, waveFunctionOld);
        } else
        {
            if (!UseImportanceSampling)
            {
                /* Brute force sampling */
                MetropolisBruteForce(rNew, rOld, waveFunctionNew, waveFunctionOld);
            } else
            {
                /* Importance sampling */
                ImportanceSampling(rNew, rOld, QForceNew, QForceOld, waveFunctionNew, waveFunctionOld);
            }
        }

        /* Write to file */
        if (cycleStepToFile != 0 && world_rank == 0 && cycle % cycleStepToFile == 0)
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
    if (cycleStepToFile != 0 && world_rank == 0)
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
        if (nParticles == 2)
        {
            waveFunctionNew = Wavefunction::TrialWaveFunction(rNew, alpha, beta, omega, spinParameter, UseJastrowFactor);
            acceptanceWeight = (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld);
        }
        else
        {
            /* the ratios acts on behalf of the wave function */
            slaterRatio = Wavefunction::SlaterRatio(rNew, nParticles, alpha, omega, InverseSlaterUpOld, InverseSlaterDownOld, i);
            jastrowRatio = Wavefunction::JastrowRatio(rNew, rOld, nParticles, beta, spinMatrix);

            acceptanceWeight = slaterRatio*slaterRatio * jastrowRatio*jastrowRatio;
        }

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
    double fraction  = 1.0/(4.0*diffusionCoefficient*timeStep);
    arma::rowvec yx  = rNew.row(i) - rOld.row(i) - diffusionCoefficient*timeStep*QForceOld.row(i);
    double yxSquared = arma::dot(yx, yx);
    arma::rowvec xy  = rOld.row(i) - rNew.row(i) - diffusionCoefficient*timeStep*QForceNew.row(i);
    double xySquared = arma::dot(xy, xy);
    return exp((-xySquared + yxSquared) * fraction) + (nParticles - 1.0);
}


void VariationalMonteCarlo::UpdateInverseSlater(const int &i)
{
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
        rOld.row(i)        = rNew.row(i);
        QForceOld.row(i)   = QForceNew.row(i);
        waveFunctionOld    = waveFunctionNew;
        UpdateInverseSlater(i);
        acceptanceCounter += 1;
    } else
    {
        rNew.row(i)      = rOld.row(i);
        QForceNew.row(i) = QForceOld.row(i);
        if (i < nParticles/2) { InverseSlaterUpNew.col(i)                = InverseSlaterUpOld.col(i); }
        else                  { InverseSlaterDownNew.col(i-nParticles/2) = InverseSlaterDownOld.col(i-nParticles/2); }
        //waveFunctionNew  = waveFunctionOld;  /* Probably unnecessary since the wavefunction didn't change. */
    }

    /* Update energies */
    if (!UseAnalyticalExpressions)
    {
        /* using numerical expressions */
        numericalEnergyVector = Hamiltonian::NumericalLocalEnergy(rNew, nParticles, nDimensions, alpha, beta, omega,
                                                                  spinParameter, UseJastrowFactor, UseFermionInteraction,
                                                                  UseNumericalPotentialEnergy);
        deltaEnergy         = numericalEnergyVector.at(0);
        kineticEnergySum   += numericalEnergyVector.at(1);
        potentialEnergySum += numericalEnergyVector.at(2);
    } else
    {
        /* using analytical expressions */
        if (nParticles == 2)
        {
            deltaEnergy = Hamiltonian::LocalEnergyTwoParticles(rNew, alpha, beta, omega, spinParameter, UseJastrowFactor,
                                                               UseFermionInteraction);
        }
        else
        {
            energyVector = Hamiltonian::LocalEnergyMoreParticles(rNew, nParticles, nDimensions, alpha, beta, omega,
                                                                spinMatrix, UseJastrowFactor, UseFermionInteraction,
                                                                InverseSlaterUpNew, InverseSlaterDownNew);
            deltaEnergy         = energyVector.at(0);
            kineticEnergySum   += energyVector.at(1);
            potentialEnergySum += energyVector.at(2);
        }

    }
    energySum        += deltaEnergy;
    energySquaredSum += (deltaEnergy*deltaEnergy);

    if (cycleType == "SteepestDescent")
    {
        if (nParticles == 2)
        {
            dPsiOfAlpha  = Wavefunction::DerivativePsiOfAlpha(rNew, omega);
            dPsiOfBeta   = Wavefunction::DerivativePsiOfBeta(rNew, beta, spinParameter);
        } else
        {
            dPsiOfAlpha  = Wavefunction::DerivativePsiManyOfAlpha(rNew, nParticles, alpha, omega);
            dPsiOfBeta   = Wavefunction::DerivativePsiManyOfBeta(rNew, nParticles, beta, spinMatrix);
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


void VariationalMonteCarlo::SteepestDescent(const int &nParticles)
{
    std::cout << "Running steepest descent..." << std::endl;
    double eta                   = 0.001;
    int nAlpha                   = 1000;
    double averagePsiAlpha       = 0.0;
    double averagePsiBeta        = 0.0;
    double averageEnergy         = 0.0;
    double alphaEnergyDerivative = 0.0;
    double betaEnergyDerivative  = 0.0;
    double scalingFactor         = 1.0/(nCycles*nParticles);
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
