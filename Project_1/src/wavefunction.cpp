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
