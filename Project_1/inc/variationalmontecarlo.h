#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>
#include <random>

class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo();
    double runMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength);
    double RandomNumber();
    void InitialTrialPositions(const arma::mat &r, int nParticles, int nDimensions);

private:
    arma::mat rOld;
    arma::mat rNew;
    arma::mat QForceOld;
    arma::mat QForceNew;
};

#endif // VARIATIONALMONTECARLO_H
