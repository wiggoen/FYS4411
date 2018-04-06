#include "inc/wavefunction.h"


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


void Wavefunction::QuantumForceInteraction(const arma::mat &r, arma::mat &QForce, double &alpha)
{
    QForce = -4.0 * alpha * r;
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
