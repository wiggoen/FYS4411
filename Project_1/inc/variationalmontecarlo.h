#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include "armadillo"


class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo();
    void runMonteCarloIntegration();

private:
    int nDimensions;
    int nParticles;
    int nCycles;
    double alpha;
    double stepLength;

    arma::mat rOld;
    arma::mat rNew;
    arma::mat QForceOld;
    arma::mat QForceNew;
};

#endif // VARIATIONALMONTECARLO_H
