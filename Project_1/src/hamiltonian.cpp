#include "inc/hamiltonian.h"
#include <armadillo>


Hamiltonian::Hamiltonian()
{

}

Hamiltonian::~Hamiltonian()
{

}

double Hamiltonian::localEnergy(const arma::mat &r, int nParticles, int nDimensions, double alpha)
{
    double localEnergy = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            localEnergy += -2 * alpha*alpha * r(i, j)*r(i, j) + alpha + 0.5 * r(i, j)*r(i, j);
        }
    }
    return localEnergy;
}
