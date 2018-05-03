#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>
#include <string>

class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo( void );
    arma::rowvec RunMonteCarloIntegration(const int nParticles, const int nDimensions, const int nCycles,
                                          const double alpha, const double stepLength,
                                          const double timeStep, const int cycleStepToFile,
                                          const std::string samplingType,
                                          const std::string derivationType,
                                          const std::string cycleType);
    void InitialTrialPositionsBruteForce(arma::mat &r);
    void InitialTrialPositionsImportanceSampling(arma::mat &r);
    void RedrawPositionImportanceSampling(arma::mat &r, const int &i);
    void CheckInitialDistance(arma::mat &rOld);
    void MonteCarloCycles( void );
    double UniformRandomNumber( void );
    double GaussianRandomNumber( void );
    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, double &waveFunctionNew,
                              const double &waveFunctionOld);
    void ImportanceSampling(arma::mat &rNew, const arma::mat &rOld, arma::mat &QForceNew,
                            const arma::mat &QForceOld, double &waveFunctionNew,
                            const double &waveFunctionOld);
    double GreensFunction(const arma::mat &rNew, const arma::mat &rOld, const arma::mat &QForceNew,
                          const arma::mat &QForceOld, const double &diffusionCoefficient,
                          const double &timeStep, const int &i);
    void UpdateEnergies(const int &i);
    double SteepestDescent(const int &nParticles);
    void SetOneBody( void );
    void OneBodyDensity( void );
    double waveFunctionOld;      /* old wave function */
    double waveFunctionNew;      /* new wave function */
    double energySum;            /* sum of particle energies for all Monte Carlo cycles */
    double energySquaredSum;     /* squared sum of particle energies for all Monte Carlo cycles */
    double deltaEnergy;          /* energy of particle we look at */
    double psiSum;
    double psiTimesEnergySum;
    double deltaPsi;
    const double a;              /* hard-core diameter of the bosons */
private:
    arma::mat rOld;              /* matrix of old position */
    arma::mat rNew;              /* matrix of new position */
    arma::mat QForceOld;         /* matrix of old quantum force */
    arma::mat QForceNew;         /* matrix of new quantum force */
    int nParticles;              /* number of particles */
    int nDimensions;             /* number of dimensions */
    int nCycles;                 /* number of Monte Carlo cycles */
    double alpha;                /* variational parameter */
    double stepLength;           /* step length used in brute force sampling */
    double timeStep;             /* time step used in importance sampling */
    int cycleStepToFile;         /* fraction of Monte Carlo cycles written to file */
    std::string samplingType;    /* sampling choice */
    std::string derivationType;  /* differentiation choice */
    std::string cycleType;       /* Monte Carlo cycles or Steepest descent */
    double beta;                 /* variational parameter */
    double acceptanceWeight;     /* weight used for accepting or rejecting a move */
    int acceptanceCounter;       /* number of accepted moves */
    arma::rowvec hist;
    arma::rowvec volume;
    double r_step;
    int nBins;
    bool isOneBodySetup;
};

#endif /* VARIATIONALMONTECARLO_H */
