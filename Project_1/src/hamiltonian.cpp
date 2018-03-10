#include "inc/hamiltonian.h"
#include <armadillo>


Hamiltonian::Hamiltonian()
{

}

Hamiltonian::~Hamiltonian()
{

}

double Hamiltonian::LocalEnergy(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha)
{
    double alphaSquared = alpha*alpha;
    double rSquared = 0;
    double localEnergy = 0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rSquared = r(i, j)*r(i, j);
            localEnergy += (-2 * alphaSquared + 0.5) * rSquared + alpha;
        }
    }
    return localEnergy;
}
