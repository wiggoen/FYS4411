#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"


Wavefunction::Wavefunction( void )
{

}


Wavefunction::~Wavefunction( void )
{

}


double Wavefunction::TrialWaveFunction(const arma::mat &r, const double &alpha, const double &beta,
                                       const double &omega, const double &a, const double &constant)
{
    double r_1Squared = r(0,0)*r(0,0) + r(0,1)*r(0,1);
    double r_2Squared = r(1,0)*r(1,0) + r(1,1)*r(1,1);
    double r_ij = abs(sqrt(r_1Squared)-sqrt(r_2Squared));
    double argument1 = -alpha*omega*(r_1Squared+r_2Squared)/2;
    double argument2 = a*r_ij/(1+beta*r_ij);
    return constant*exp(argument1+argument2);
}
