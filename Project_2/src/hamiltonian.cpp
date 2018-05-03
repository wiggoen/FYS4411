#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include <armadillo>


Hamiltonian::Hamiltonian( void )
{

}


Hamiltonian::~Hamiltonian( void )
{

}


double Hamiltonian::LocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions, const double &alpha)
{

    return localEnergy;
}

