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
    arma::rowvec RunVMC(const int nParticles, const int nCycles, const double alpha,
                        const double beta,  const double omega, const double spinParameter,
                        const double stepLength, const double timeStep, const bool UseJastrowFactor,
                        const bool UseImportanceSampling, const bool UseFermionInteraction,
                        std::string derivationType, std::string cycleType);
    void InitialTrialPositionsBruteForce(arma::mat &r);
    void MonteCarloCycles( void );
    double UniformRandomNumber( void );
    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, double &waveFunctionNew,
                                 const double &waveFunctionOld);
    void ImportanceSampling(arma::mat &rNew, const arma::mat &rOld, arma::mat &QForceNew,
                                 const arma::mat &QForceOld, double &waveFunctionNew,
                                 const double &waveFunctionOld);
    void UpdateEnergies(const int &i);
    void InitialTrialPositionsImportanceSampling(arma::mat &r);
    double GaussianRandomNumber( void );
    double GreensFunction(const arma::mat &rNew, const arma::mat &rOld, const arma::mat &QForceNew,
                          const arma::mat &QForceOld, const double &diffusionCoefficient,
                          const double &timeStep, const int &i);
    double waveFunctionOld;      /* old wave function */
    double waveFunctionNew;      /* new wave function */
    double energySum;            /* sum of particle energies for all Monte Carlo cycles */
    double energySquaredSum;     /* squared sum of particle energies for all Monte Carlo cycles */
    double deltaEnergy;          /* energy of particle we look at */
    //double psiSum;
    //double psiTimesEnergySum;
    //double deltaPsi;
    const int nDimensions;       /* number of dimensions */
private:
    arma::mat rOld;              /* matrix of old position */
    arma::mat rNew;              /* matrix of new position */
    arma::mat QForceOld;         /* matrix of old quantum force */
    arma::mat QForceNew;         /* matrix of new quantum force */
    int nParticles;              /* number of particles */
    int nCycles;                 /* number of Monte Carlo cycles */
    double alpha;                /* variational parameter */
    double beta;                 /* variational parameter */
    double omega;
    double spinParameter;        /* spin parameter a */
    double stepLength;           /* step length used in brute force sampling */
    double timeStep;             /* time step used in importance sampling */
    bool UseJastrowFactor;       /* Wave function with/without Jastrow factor */
    bool UseImportanceSampling;  /* With/without importance sampling */
    bool UseFermionInteraction;  /* With/without interaction */
    std::string derivationType;  /* differentiation choice */
    std::string cycleType;       /* Monte Carlo cycles or Steepest descent */
    //int cycleStepToFile;       /* fraction of Monte Carlo cycles written to file */
    //std::string samplingType;    /* sampling choice */
    double acceptanceWeight;     /* weight used for accepting or rejecting a move */
    int acceptanceCounter;       /* number of accepted moves */
    //arma::rowvec hist;
    //arma::rowvec volume;
    //double r_step;
    //int nBins;
    //bool isOneBodySetup;

};

#endif /* VARIATIONALMONTECARLO_H */
