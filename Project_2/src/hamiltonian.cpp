#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include <armadillo>


Hamiltonian::Hamiltonian( void )
{

}


Hamiltonian::~Hamiltonian( void )
{

}


double Hamiltonian::LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha,
                                const double &beta, const double &omega, const double &a)
{
    double localEnergy = 0;
    double twoAlphaOmega = alpha*omega;         // REMOVE two??!?!
    double omegaSquaredHalf = 0.5*omega*omega;

    double betaR_ij =  beta*arma::norm(r.row(0) - r.row(1));
    //std::cout << "norm = " << betaR_ij << std::endl;

    for (int i = 0; i < nParticles; i++)
    {
        double numerator = a*(1 + betaR_ij) - a*betaR_ij;
        double denominator = (1 + betaR_ij)*(1 + betaR_ij);
        double firstFraction = numerator/denominator;
        //double firstTerm = -0.5*(-twoAlphaOmega*arma::norm(r.row(i)) + firstFraction);
        double firstTerm = -twoAlphaOmega*arma::norm(r.row(i)) + firstFraction;

        //double secondTerm = 1 - twoAlphaOmega*arma::dot(r.row(i), r.row(i)) + firstFraction*arma::norm(r.row(i));

        double secondTerm = firstFraction * (2*beta*arma::norm(r.row(i)) + 2*beta*beta*arma::norm(r.row(i)))/denominator - twoAlphaOmega;

        //double thirdTerm = arma::norm(r.row(i))*(-twoAlphaOmega + (((1 + beta - a*beta)*(1 + betaR_ij) - numerator)*2*beta)/(denominator*(1 + betaR_ij)));

        double thirdTerm = omegaSquaredHalf*arma::dot(r.row(i), r.row(i));

        //double fourthTerm = omegaSquaredHalf*arma::dot(r.row(i), r.row(i));

        localEnergy += -0.5*firstTerm*firstTerm * (-0.5)*secondTerm * thirdTerm;
    }

    return localEnergy;
}

