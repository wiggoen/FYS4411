#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>
#include <string>

class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo( void );
    arma::rowvec RunMonteCarloIntegration(const int nParticles, const int nCycles, const double alpha,
                                          const double beta,  const double omega, const double a, const double stepLength);
    void InitialTrialPositionsBruteForce(arma::mat &r);
    void MonteCarloCycles( void );
    double UniformRandomNumber( void );
    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, double &waveFunctionNew,
                              const double &waveFunctionOld);
    void UpdateEnergies(const int &i);
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
    //arma::mat QForceOld;         /* matrix of old quantum force */
    //arma::mat QForceNew;         /* matrix of new quantum force */
    int nParticles;              /* number of particles */
    int nCycles;                 /* number of Monte Carlo cycles */
    double alpha;                /* variational parameter */
    double stepLength;           /* step length used in brute force sampling */
    double constant;             /* normalization factor */
    //double timeStep;             /* time step used in importance sampling */
    //int cycleStepToFile;         /* fraction of Monte Carlo cycles written to file */
    std::string samplingType;    /* sampling choice */
    std::string derivationType;  /* differentiation choice */
    std::string cycleType;       /* Monte Carlo cycles or Steepest descent */
    double beta;                 /* variational parameter */
    double acceptanceWeight;     /* weight used for accepting or rejecting a move */
    int acceptanceCounter;       /* number of accepted moves */
    double a;              /* spinn parameter */
    double omega;
    //arma::rowvec hist;
    //arma::rowvec volume;
    //double r_step;
    //int nBins;
    //bool isOneBodySetup;

};

#endif /* VARIATIONALMONTECARLO_H */
