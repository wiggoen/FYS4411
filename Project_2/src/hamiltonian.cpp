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

double Hamiltonian::NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                         const double &alpha, const double &stepLength, const double &beta,
                                         const double &omega, const double &a, const double &constant)
{
    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0.0;
    double waveFunctionPlus  = 0.0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, alpha, beta, omega, a, constant);

    double stepLengthSquaredFraction = 1.0 / (stepLength * stepLength);

    /* Kinetic energy */
    double kineticEnergy = 0.0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i, j) += stepLength;
            rMinus(i, j) -= stepLength;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, alpha, beta, omega, a, constant);
            waveFunctionPlus = Wavefunction::TrialWaveFunction(rPlus, alpha, beta, omega, a, constant);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2.0 * waveFunctionCurrent);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    kineticEnergy = 0.5 * stepLengthSquaredFraction * kineticEnergy / waveFunctionCurrent;

    /* External potential */
    double externalPotential = 0.0;
    double rSquared = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        rSquared = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            rSquared += r(i, j) * r(i, j);
        }
        externalPotential += 0.5 * rSquared;
    }
    return kineticEnergy + externalPotential;
}
