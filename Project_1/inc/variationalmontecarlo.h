#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include <armadillo>
#include <string>

class VariationalMonteCarlo
{
public:
    VariationalMonteCarlo();
    ~VariationalMonteCarlo();
    arma::rowvec RunMonteCarloIntegration(const int, const int, const int, const double, const double, const double, const int, const std::string, const std::string, const std::string);
    void InitialTrialPositionsBruteForce(arma::mat &);
    void InitialTrialPositionsImportanceSampling(arma::mat &);
    void RedrawPositionImportanceSampling(arma::mat &, const int &);
    void CheckInitialDistance(arma::mat &);
    void MonteCarloCycles();
    double UniformRandomNumber();
    double GaussianRandomNumber();
    void MetropolisBruteForce(arma::mat &, arma::mat &, double &, const double &);
    void ImportanceSampling(arma::mat &, const arma::mat &, arma::mat &, const arma::mat &, double &, const double &);
    double GreensFunction(const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const double &, const double &, const int &);
    void UpdateEnergies(const int &);
    double SteepestDescent(const int &);
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
    std::string integrationType;
    std::string cycleType;
    double beta;
    double acceptanceWeight;
    int acceptanceCounter;
    /*
    std::string strParticles;			// <<< REMOVE THESE!!
    std::string strDimensions;
    std::string strAlpha;
    */
};

#endif // VARIATIONALMONTECARLO_H
