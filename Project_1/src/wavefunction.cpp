#include "inc/wavefunction.h"


Wavefunction::Wavefunction()
{

}


Wavefunction::~Wavefunction()
{

}


double Wavefunction::TrialWaveFunction(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha)
{
    double argument = 0;
    for (int i = 0; i < nParticles; i++)
    {
        double rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++)
        {
            rSingleParticle += r(i, j) * r(i, j);
        }
        argument += rSingleParticle;
    }
    return exp(-argument * alpha);
}


void Wavefunction::QuantumForce(const arma::mat &r, arma::mat &QForce, double &alpha)
{
    QForce = -4.0 * alpha * r;
}

double Wavefunction::derivativePsi(const arma::mat &R, int nParticles, int nDimensions, double beta)
// Returns 1/psi * psi'
{
    double derivative = 0;
    for (int i = 0; i<nParticles; i++)
    {
        for (int j=0; j<nDimensions; j++)
        {
            if (j<2) {derivative -=      R(i,j)*R(i,j);}
            else     {derivative -= beta*R(i,j)*R(i,j);}
        }
    }
    return derivative;
}
