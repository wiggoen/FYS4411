#pragma once
#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>
#include <string>

class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo( void );
    arma::rowvec RunVMC(const int nParticles, const int nCycles, const double alpha, const double beta, const double omega,
                        const double stepLength, const double timeStep, const bool UseJastrowFactor,
                        const bool UseImportanceSampling, const bool UseFermionInteraction,
                        const bool UseAnalyticalExpressions, const bool UseNumericalPotentialEnergy,
                        const std::string cycleType, const int cycleStepToFile, int terminalizationFactor);

    void InitialTrialPositionsBruteForce(arma::mat &r);

    void InitialTrialPositionsImportanceSampling(arma::mat &r);

    void SlaterInitialization( void );

    double UniformRandomNumber( void );

    double GaussianRandomNumber( void );

    void MonteCarloCycles( void );

    int PickRandomParticle( void );

    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld);

    void ImportanceSampling(arma::mat &rNew, const arma::mat &rOld, arma::mat &QForceNew, const arma::mat &QForceOld);

    double GreensRatio(const arma::mat &rNew, const arma::mat &rOld, const arma::mat &QForceNew, const arma::mat &QForceOld,
                       const double &diffusionCoefficient, const double &timeStep, const int &i);

    void UpdateInverseSlater(const int &i);

    void UpdateEnergies(const int &i);

    void SteepestDescent( void );

    void SetOneBody( void );

    void OneBodyDensity( void );

    double waveFunctionOld;            /* old wave function */
    double waveFunctionNew;            /* new wave function */
    double energySum;                  /* sum of particle energies for all Monte Carlo cycles */
    double energySquaredSum;           /* squared sum of particle energies for all Monte Carlo cycles */
    double deltaEnergy;                /* energy of particle we look at */
    arma::rowvec energyVector;
    double kineticEnergySum;
    double potentialEnergySum;
    double kineticEnergy;
    double potentialEnergy;
    double psiSum;
    double psiTimesEnergySum;
    const int nDimensions;             /* number of dimensions */
    double dPsiOfAlpha;
    double dPsiOfBeta;
    double psiSumAlpha;
    double psiSumBeta;
    double psiTimesEnergyA;
    double psiTimesEnergyB;
    double psiOfAlphaTimesEnergySum;
    double psiOfBetaTimesEnergySum;
    int world_rank;                    /* to use with MPI; holds the processor rank */
    arma::mat SlaterUpOld;
    arma::mat SlaterUpNew;
    arma::mat SlaterDownOld;
    arma::mat SlaterDownNew;
    arma::mat InverseSlaterUpOld;
    arma::mat InverseSlaterUpNew;
    arma::mat InverseSlaterDownOld;
    arma::mat InverseSlaterDownNew;

private:
    arma::mat rOld;                    /* matrix of old position */
    arma::mat rNew;                    /* matrix of new position */
    arma::mat QForceOld;               /* matrix of old quantum force */
    arma::mat QForceNew;               /* matrix of new quantum force */
    int nParticles;                    /* number of particles */
    int nCycles;                       /* number of Monte Carlo cycles */
    double alpha;                      /* variational parameter */
    double beta;                       /* variational parameter */
    double omega;
    double stepLength;                 /* step length used in brute force sampling */
    double timeStep;                   /* time step used in importance sampling */
    bool UseJastrowFactor;             /* Wave function with/without Jastrow factor */
    bool UseImportanceSampling;        /* With/without importance sampling */
    bool UseFermionInteraction;        /* With/without interaction */
    bool UseAnalyticalExpressions;     /* differentiation choice */
    bool UseNumericalPotentialEnergy;  /* Choice between numerical kinetic energy or full numerical energy */
    std::string cycleType;             /* Monte Carlo cycles or Steepest descent */
    int cycleStepToFile;               /* fraction of Monte Carlo cycles written to file */
    double acceptanceWeight;           /* weight used for accepting or rejecting a move */
    int acceptanceCounter;             /* number of accepted moves */
    int cycleNumber;
    int terminalizationFactor;         /* throws away the first (int) energies */
    double slaterRatio;
    double jastrowRatio;
    arma::rowvec hist;
    arma::rowvec volume;
    double r_step;
    int nBins;
    bool isOneBodySetup;
};

#endif /* VARIATIONALMONTECARLO_H */
