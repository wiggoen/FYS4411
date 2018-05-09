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

void Wavefunction::NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                         const int &nDimensions, const double &alpha, const double &stepLength,
                                         const double &beta,const double &omega, const double &a, const double &constant)
{
    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0.0;
    double waveFunctionPlus = 0.0;
    double waveFunctionCurrent = TrialWaveFunction(r, alpha, beta, omega, a, constant);

    double stepLengthFraction = 1.0 / stepLength;

    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i, j) += stepLength;
            rMinus(i, j) -= stepLength;
            waveFunctionMinus = TrialWaveFunction(rMinus, alpha, beta, omega, a, constant);
            waveFunctionPlus = TrialWaveFunction(rPlus, alpha, beta, omega, a, constant);
            QForce(i, j) = (waveFunctionPlus - waveFunctionMinus);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    QForce = QForce * stepLengthFraction / waveFunctionCurrent;
}
