#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>


class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo();
    arma::rowvec RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength, double dt, int cycleStepToFile);
    void InitialTrialPositions(arma::mat &r);
    void MonteCarloCycles();
    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew);
    double UniformRandomNumber();
    double GaussianRandomNumber();
    void FokkerPlanckAndLangevin(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew);
    double GreensFunction(double &oldPosition, double &newPosition, double &D, double &dt, double &QForceOld);
    void UpdateEnergies();
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
    double dt;
    int cycleStepToFile;
    int x;
    int y;
    double acceptanceWeight;
    double acceptanceCounter;
};

#endif // VARIATIONALMONTECARLO_H
