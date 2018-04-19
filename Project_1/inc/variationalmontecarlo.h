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
    double waveFunctionOld;
    double waveFunctionNew;
    double energySum;
    double energySquaredSum;
    double deltaEnergy;
    double psiSum;
    double psiTimesEnergySum;
    double deltaPsi;
    const double a;
private:
    arma::mat rOld;
    arma::mat rNew;
    arma::mat QForceOld;
    arma::mat QForceNew;
    int nParticles;
    int nDimensions;
    int nCycles;
    double alpha;
    double stepLength;
    double timeStep;
    int cycleStepToFile;
    std::string samplingType;
    std::string derivationType;
    std::string cycleType;
    double beta;
    double acceptanceWeight;
    int acceptanceCounter;
    int thrownCounter;
};

#endif /* VARIATIONALMONTECARLO_H */
