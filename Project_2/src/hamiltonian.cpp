#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include <armadillo>


Hamiltonian::Hamiltonian( void )
{

}


Hamiltonian::~Hamiltonian( void )
{

}


double Hamiltonian::LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha, const double &beta, const double &omega, const double &a)
{
    double localEnergy = 0;
    double twoAlphaOmega = 2*alpha*omega;
    double omegaSquaredHalf = 0.5*omega*omega;

    //for (int i = 0; i < nParticles; i++)
    int i = 0;
    {
        double betaR_ij =  beta*arma::norm(r.row(i) - r.row(i+1));
        double numerator = 1 + betaR_ij - a*betaR_ij;
        double denominator = (1 + betaR_ij)*(1 + betaR_ij);
        double firstFraction = numerator/denominator;
        double firstTerm = -0.5*(-twoAlphaOmega*arma::norm(r.row(i)) + firstFraction);

        double secondTerm = 1 - twoAlphaOmega*arma::dot(r.row(i), r.row(i)) + firstFraction*arma::norm(r.row(i));

        double thirdTerm = arma::norm(r.row(i))*(-twoAlphaOmega + (((1 + beta - a*beta)*(1 + betaR_ij) - numerator)*2*beta)/(denominator*(1 + betaR_ij)));

        double fourthTerm = omegaSquaredHalf*arma::dot(r.row(i), r.row(i));

        localEnergy += firstTerm * secondTerm * thirdTerm * fourthTerm;
    }

    return localEnergy;
}

