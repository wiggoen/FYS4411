#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>
#include <random>

class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo();
    double RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength, int nDataPoints);
    double RandomNumber();
    double GaussianRandomNumber();
    void InitialTrialPositions(arma::mat &r);
    void MonteCarloCycles();
    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew, int i);
    void FokkerPlanckAndLangevin(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew, int i);
    double GreensFunction(double x, double y, double D, double deltaT, double QForceOld);
    double waveFunctionOld;
    double waveFunctionNew;
    double energySum;
    double energySquaredSum;
    double deltaEnergy;
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
    int nDataPoints;
};

#endif // VARIATIONALMONTECARLO_H
