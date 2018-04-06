#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include <armadillo>
#include <limits>               // trenger vi denne? infinity

Hamiltonian::Hamiltonian()
{

}


Hamiltonian::~Hamiltonian()
{

}


double Hamiltonian::LocalEnergy(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha)
{
    double factors = -2.0 * alpha*alpha + 0.5;
    double dimensionality = nDimensions * alpha;

    double localEnergy = 0;
    for (int i = 0; i < nParticles; i++)
    {
        double rSquared = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            rSquared += r(i, j)*r(i, j);
        }
        localEnergy += factors * rSquared + dimensionality;
    }
    return localEnergy;
}


double Hamiltonian::NumericalLocalEnergy(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha, double &stepLength)
{
    double beta = 1;

    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, nParticles, nDimensions, alpha, beta);

    double stepLengthSquaredFraction = 1.0 / (stepLength * stepLength);

    // Kinetic energy
    double kineticEnergy = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i, j) += stepLength;
            rMinus(i, j) -= stepLength;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, nParticles, nDimensions, alpha, beta);
            waveFunctionPlus = Wavefunction::TrialWaveFunction(rPlus, nParticles, nDimensions, alpha, beta);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    kineticEnergy = 0.5 * stepLengthSquaredFraction * kineticEnergy / waveFunctionCurrent;

    // External potential
    double externalPotential = 0;
    double rSquared = 0;
    for (int i = 0; i < nParticles; i++)
    {
        rSquared = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            rSquared += r(i, j) * r(i, j);
        }
        externalPotential += 0.5 * rSquared;
    }
    return kineticEnergy + externalPotential;
}


double Hamiltonian::LocalEnergyInteraction(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha, double &beta)
{
    //std::cout << "beta = " << beta << std::endl;

    double a = 0.0043;
    arma::mat PhiGradient = -2 * alpha * r;
    double PhiLaplacian = 0;
    arma::rowvec vectorSum;
    double gradientProduct = 0;
    double doubleSum = 0;
    double derivativeSum = 0;
    double singleParticle = 0;
    //double distanceFractionSum = 0;

    double repulsivePontential = 0;

    double laplacianSum = 0;
    for (int k = 0; k < nParticles; k++)
    {
        PhiLaplacian = arma::dot(PhiGradient.row(k), PhiGradient.row(k)) - 2 * alpha * (2 + beta);

        vectorSum = VectorSum(r, nParticles, nDimensions, a, k);
        gradientProduct = arma::dot(PhiGradient.row(k), vectorSum);
        doubleSum = arma::dot(vectorSum, vectorSum);
        derivativeSum = DerivativeSum(r, nParticles, a, k);

        //std::cout << "gradientProduct = " << gradientProduct << std::endl;
        //std::cout << "doubleSum = " << doubleSum << std::endl;

        singleParticle += arma::dot(r.row(k), r.row(k));
        //distanceFractionSum += Correlation(r, nParticles, a, k);          << Trengs denne?

        repulsivePontential += RepulsivePotential(r, nParticles, a, k);     // << er dette korrekt implementasjon?

        laplacianSum -= (PhiLaplacian + gradientProduct + doubleSum + derivativeSum);
    }


    // TODO: Implement!!
    // Repulsive potential
    /*
    double rij = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            rij = 0;
            for(int k = 0; k < nDimensions; k++) {
                rij += (r(i, k) - r(j, k)) * (r(i, k) - r(j, k));
            }
            potentialEnergy += 1 / sqrt(rij);
        }
    }
    */
    std::cout << "laplacianSum = " << laplacianSum << std::endl;
    return 0.5 * (laplacianSum + singleParticle * exp(-singleParticle * alpha)) + repulsivePontential;// * distanceFractionSum;//0.5 * (f+g+h) + pot;
}


double Hamiltonian::ParticleDistance(const arma::rowvec &r_k, const arma::rowvec &r_j)
{
    double distance = arma::norm(r_k - r_j);
    return distance;
}


arma::rowvec Hamiltonian::VectorSum(const arma::mat &r, int &nParticles, int &nDimensions, const double a, const int &k)
{
    arma::rowvec vectorSum = arma::zeros<arma::rowvec>(nDimensions);
    double distance = 0;

    for (int j = 0; j < nParticles; j++)
    {
        if (j != k)
        {
            //std::cout << "j = " << j << "     " << "k = " << k << std::endl;
            distance = ParticleDistance(r.row(k), r.row(j));
            vectorSum += (((r.row(k) - r.row(j)) / distance) * (a / ((distance*distance)*(distance-a))));
        }
    }
    //std::cout << "vectorSum = " << vectorSum << std::endl;
    return vectorSum;
}


double Hamiltonian::DerivativeSum(const arma::mat &r, int &nParticles, const double a, const int &k)
{
    double derivativeSum = 0;
    double distance = 0;

    for (int j = 0; j < nParticles; j++)
    {
        if (j != k)
        {
            //std::cout << "j = " << j << "     " << "k = " << k << std::endl;
            distance = ParticleDistance(r.row(k), r.row(j));
            derivativeSum += (-a*a / (distance*distance * (distance - a)*(distance - a)));
        }
    }
    //std::cout << "derivativeSum = " << derivativeSum << std::endl;
    return derivativeSum;
}


double Hamiltonian::Correlation(const arma::mat &r, int &nParticles, double &a, int &k)
{
    double distanceFractionSum = 0;

    for (int j = 0; j < nParticles; j++)
    {
        if (j != k)
        {
            distanceFractionSum += (1 - a/ParticleDistance(r.row(k), r.row(j)));
        }
    }
    return distanceFractionSum;
}


double Hamiltonian::RepulsivePotential(const arma::mat &r, int &nParticles, double &a, int &k)
{
    double repulsivePotential = 0;
    double distance = 0;
    double infinity = std::numeric_limits<double>::infinity();

    //std::cout << "infinity = " << infinity << std::endl;

    for (int j = 0; j < nParticles; j++)
    {
        if (k < j)
        {
            distance = ParticleDistance(r.row(k), r.row(j));
            if (distance <= a)
            {
                repulsivePotential += infinity;      // UENDELIG?!?!?
            } else
            {
                repulsivePotential += 0;
            }
        }
    }
    return repulsivePotential;
}
