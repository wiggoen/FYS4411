#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"


Wavefunction::Wavefunction()
{

}


Wavefunction::~Wavefunction()
{

}


double Wavefunction::TrialWaveFunction(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha, double &beta)
{
    double argument = 0;
    for (int i = 0; i < nParticles; i++)
    {
        double rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            if (j < 2)
            {
                rSingleParticle += r(i, j) * r(i, j);
            } else
            {
                rSingleParticle += beta * r(i, j) * r(i, j);
            }
        }
        argument += rSingleParticle;
    }
    return exp(-argument * alpha);
}


double Wavefunction::TrialWaveFunctionInteraction(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha, double &beta, double &a)
{
    double trailWaveFunctionOneBody = Wavefunction::TrialWaveFunction(r, nParticles, nDimensions, alpha, beta);
    double distance = 0;
    double f = 1;

    for (int i = 0; i < nParticles; i++)
    {
        for (int j = i+1; j < nParticles; j++)
        {
            distance = Hamiltonian::ParticleDistance(r.row(i), r.row(j));

            if (distance > a) {
                f *= (1 - a / distance);
            } else
            {
                f *= 0;
            }
        }
    }
    return trailWaveFunctionOneBody * f;
}


void Wavefunction::QuantumForce(const arma::mat &r, arma::mat &QForce, double &alpha)
{
    QForce = -4.0 * alpha * r;
}


void Wavefunction::NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, int &nParticles, int &nDimensions, double &alpha, double &stepLength)
{
    double beta = 1;

    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;
    double waveFunctionCurrent = TrialWaveFunction(r, nParticles, nDimensions, alpha, beta);

    double stepLengthFraction = 1.0 / stepLength;

    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i, j) += stepLength;
            rMinus(i, j) -= stepLength;
            waveFunctionMinus = TrialWaveFunction(rMinus, nParticles, nDimensions, alpha, beta);
            waveFunctionPlus = TrialWaveFunction(rPlus, nParticles, nDimensions, alpha, beta);
            QForce(i, j) = (waveFunctionPlus - waveFunctionMinus);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    QForce = QForce * stepLengthFraction / waveFunctionCurrent;
}


void Wavefunction::QuantumForceInteraction(const arma::mat &r, arma::mat &QForce, double &alpha, int &nParticles, int &nDimensions, double &a, int &i)
{
    arma::rowvec vectorSum = Hamiltonian::VectorSum(r, nParticles, nDimensions, a, i);
    QForce.row(i) = -4.0 * alpha * r.row(i) + 2 * a * vectorSum;
}


// TODO: CHECK IF THIS IS IN USE.
double Wavefunction::DerivativePsi(const arma::mat &r, int &nParticles, int &nDimensions, double &beta)
// Returns 1/psi * psi'
{
    double derivative = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            if (j < 2) {
                derivative -= r(i,j) * r(i,j);
            }
            else
            {
                derivative -= beta * r(i,j) * r(i,j);
            }
        }
    }
    return derivative;
}
