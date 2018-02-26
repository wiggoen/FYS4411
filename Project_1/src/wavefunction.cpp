#include "inc/wavefunction.h"
#include "inc/matrix.h"
#include <math.h>
#include <cmath>

Wavefunction::Wavefunction(const arma::mat &r)
{

}

double Wavefunction::waveFunction(const arma::mat &r)
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

void Wavefunction::QuantumForce(const arma::mat &r, arma::mat &QForce)
{
    QForce = -4 * alpha * r;
}
