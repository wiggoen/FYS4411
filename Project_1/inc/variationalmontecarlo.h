#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>
#include <string>

class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo();
    arma::rowvec RunMonteCarloIntegration(int nParticles, int nDimensions, int nCycles, double alpha, double stepLength, double timeStep, int cycleStepToFile);
    void InitialTrialPositionsBruteForce(arma::mat &r);
    void InitialTrialPositionsImportanceSampling(arma::mat &r);
    void MonteCarloCycles();
    double UniformRandomNumber();
    double GaussianRandomNumber();
    void MetropolisBruteForce(arma::mat &rNew, arma::mat &rOld, double &waveFunctionOld, double &waveFunctionNew);
    void ImportanceSampling(arma::mat &rNew, arma::mat &rOld, arma::mat &QForceOld, arma::mat &QForceNew, double &waveFunctionOld, double &waveFunctionNew);
    double GreensFunction(const arma::mat &rOld, const arma::mat &rNew, const arma::mat &QForceOld, double &diffusionCoefficient, double &timeStep, int &i);
    void UpdateEnergies(int &i);
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
    double timeStep;
    int cycleStepToFile;
    double acceptanceWeight;
    int acceptanceCounter;
    std::string samplingType;
	/*
	std::string strParticles;			// <<< REMOVE THESE!!
    std::string strDimensions;
    std::string strAlpha;
	*/
};

#endif // VARIATIONALMONTECARLO_H
