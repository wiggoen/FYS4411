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
