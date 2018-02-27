#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>
#include <random>

class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo();
    double RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength);
    double RandomNumber();
    void InitialTrialPositions(arma::mat &r);
    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew, int i);
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
};

#endif // VARIATIONALMONTECARLO_H
