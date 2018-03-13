#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>


class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo();
    arma::rowvec RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength, double dt, int cycleStepToFile);
    void InitialTrialPositionsBruteForce(arma::mat &r);
    void InitialTrialPositionsImportanceSampling(arma::mat &r);
    void MonteCarloCycles();
    double UniformRandomNumber();
    double GaussianRandomNumber();
    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, double &waveFunctionOld, double &waveFunctionNew);
    void ImportanceSampling(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew);
    double GreensFunction(const arma::mat &rOld, const arma::mat &rNew, const arma::mat &QForceOld, double &D, double &dt, int &i);
    void UpdateEnergies(int &i);
    void Blocking(arma::vec samples, int nSamples, int block_size, arma::vec results);
    double Mean(arma::vec &samples, int nSamples);
    double Variance(arma::vec &samlpes, int nSamples, arma::vec &results);
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
    double acceptanceWeight;
    double acceptanceCounter;
};

#endif // VARIATIONALMONTECARLO_H
