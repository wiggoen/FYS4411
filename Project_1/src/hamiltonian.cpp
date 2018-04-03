#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include <armadillo>


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
    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, nParticles, nDimensions, alpha);

    double stepLengthSquaredFraction = 1.0 / (stepLength * stepLength);

    // Kinetic energy
    double kineticEnergy = 0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i, j) += stepLength;
            rMinus(i, j) -= stepLength;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, nParticles, nDimensions, alpha);
            waveFunctionPlus = Wavefunction::TrialWaveFunction(rPlus, nParticles, nDimensions, alpha);
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

    // TODO: Implement!!
    // Internal potential
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
    return kineticEnergy + externalPotential;
}
