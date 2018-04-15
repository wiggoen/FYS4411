#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include <armadillo>


Hamiltonian::Hamiltonian()
{

}


Hamiltonian::~Hamiltonian()
{

}


double Hamiltonian::LocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions, const double &alpha)
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


double Hamiltonian::NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                         const double &alpha, const double &stepLength, const double &beta)
{
    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0.0;
    double waveFunctionPlus  = 0.0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, nParticles, nDimensions, alpha, beta);

    double stepLengthSquaredFraction = 1.0 / (stepLength * stepLength);

    // Kinetic energy
    double kineticEnergy = 0.0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i, j) += stepLength;
            rMinus(i, j) -= stepLength;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, nParticles, nDimensions, alpha, beta);
            waveFunctionPlus = Wavefunction::TrialWaveFunction(rPlus, nParticles, nDimensions, alpha, beta);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2.0 * waveFunctionCurrent);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    kineticEnergy = 0.5 * stepLengthSquaredFraction * kineticEnergy / waveFunctionCurrent;

    // External potential
    double externalPotential = 0.0;
    double rSquared = 0.0;
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


double Hamiltonian::LocalEnergyInteraction(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                           const double &alpha, const double &beta, const double &a)
{
    double aSquared = a*a;
    double betaSquared = beta*beta;

    double laplacianTerm = 0.0;                // first term of the second derivative of the wavefunction
    double gradientTerm  = 0.0;                // second term of the second derivative of the wavefunction
    double doubleSum     = 0.0;                // third term of the second derivative of the wavefunction
    double derivativeSum = 0.0;                // fourth term of the second derivative of the wavefunction
    double secondDerivativeOfWavefunction = 0.0;

    arma::mat PhiGradient = -2.0 * alpha * r;  // the first part of the gradient term
    if (nDimensions > 2)
    {
        for (int i = 0; i < nParticles; i++)
        {
            PhiGradient(i, 2) *= beta;
        }
    }
    arma::rowvec vectorSum;                    // the second part of the gradient term (sum in the gradient term)

    double externalPotential  = 0.0;
    double repulsivePotential = 0.0;

    for (int k = 0; k < nParticles; k++)
    {
        laplacianTerm = arma::dot(PhiGradient.row(k), PhiGradient.row(k)) - 2.0 * alpha * (2.0 + beta);

        vectorSum = VectorSum(r, nParticles, nDimensions, a, k);
        gradientTerm = 2.0 * a * arma::dot(PhiGradient.row(k), vectorSum);

        derivativeSum = - aSquared * DerivativeSum(r, nParticles, a, k);

        doubleSum = derivativeSum*derivativeSum;

        secondDerivativeOfWavefunction += (laplacianTerm + gradientTerm + doubleSum + derivativeSum);

        repulsivePotential += RepulsivePotential(r, nParticles, a, k);

        double argument = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            if (j < 2) { argument += r(k, j) * r(k, j); }
            else       { argument += betaSquared * r(k, j) * r(k, j); }
        }
        externalPotential += argument;
    }
    return 0.5 * (-secondDerivativeOfWavefunction + externalPotential) + repulsivePotential;
}


double Hamiltonian::ParticleDistance(const arma::rowvec &r_k, const arma::rowvec &r_j)
{
    return arma::norm(r_k - r_j);
}


arma::rowvec Hamiltonian::VectorSum(const arma::mat &r, const int &nParticles, const int &nDimensions, const double &a,
                                    const int &k)
{
    arma::rowvec vectorSum = arma::zeros<arma::rowvec>(nDimensions);
    double distance = 0.0;
    double distanceSquared = 0.0;

    for (int j = 0; j < nParticles; j++)
    {
        if (j != k)
        {
            distance = ParticleDistance(r.row(k), r.row(j));
            distanceSquared = distance*distance;
            vectorSum += ( (r.row(k) - r.row(j)) / (distanceSquared*(distance - a)) );
        }
    }
    return vectorSum;
}


double Hamiltonian::DerivativeSum(const arma::mat &r, const int &nParticles, const double &a, const int &k)
{
    double derivativeSum = 0.0;
    double distance = 0.0;
    double distanceSquared = 0.0;
    double distanceSubtraction = 0.0;          // distance - a
    double distanceSubtractionSquared = 0.0;

    for (int j = 0; j < nParticles; j++)
    {
        if (j != k)
        {
            distance = ParticleDistance(r.row(k), r.row(j));
            distanceSquared = distance*distance;
            distanceSubtraction = distance - a;
            distanceSubtractionSquared = distanceSubtraction*distanceSubtraction;
            derivativeSum += (1.0 / (distanceSquared * distanceSubtractionSquared) );
        }
    }
    return derivativeSum;
}


double Hamiltonian::RepulsivePotential(const arma::mat &r, const int &nParticles, const double &a, const int &k)
{
    double repulsivePotential = 0.0;
    double distance = 0.0;

    for (int j = 0; j < nParticles; j++)
    {
        if (k < j)
        {
            distance = ParticleDistance(r.row(k), r.row(j));
            if (distance <= a)
            {
                repulsivePotential += 1e16;
            }
        }
    }
    return repulsivePotential;
}
