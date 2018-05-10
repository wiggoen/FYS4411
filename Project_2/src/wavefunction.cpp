#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"


Wavefunction::Wavefunction( void )
{

}


Wavefunction::~Wavefunction( void )
{

}


double Wavefunction::TrialWaveFunction(const arma::mat &r, const double &alpha, const double &beta,
                                       const double &omega, const double &spinParameter, const bool &UseJastrowFactor)
{
    double r_1Squared = r(0,0)*r(0,0) + r(0,1)*r(0,1);
    double r_2Squared = r(1,0)*r(1,0) + r(1,1)*r(1,1);
    double argument1 = -0.5*alpha*omega*(r_1Squared + r_2Squared);

    if (!UseJastrowFactor) {
        return exp(argument1);
    } else {
        double r_ij = fabs(sqrt(r_1Squared)-sqrt(r_2Squared));
        double argument2 = spinParameter*r_ij/(1+beta*r_ij);
        return exp(argument1+argument2);
    }
}

void Wavefunction::QuantumForce(const arma::mat &r, arma::mat &QForce, const double &alpha,
                                const double &beta, const double &omega, const double &spinParameter)
{
    double x1 = r(0,0); double x2 = r(1,0);
    double y1 = r(0,1); double y2 = r(1,1);
    double r12 = abs(sqrt(r(0,0)*r(0,0)+r(0,1)*r(0,1))-sqrt(r(1,0)*r(1,0)+r(1,1)*r(1,1)));
    QForce =  -alpha*omega*r - spinParameter*(x1-x2)/(r12*(1+beta*r12)*(1+beta*r12));
    QForce += -alpha*omega*r - spinParameter*(y1-y2)/(r12*(1+beta*r12)*(1+beta*r12));
}


/*
void Wavefunction::NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                         const int &nDimensions, const double &alpha, const double &beta, const double &omega,
                                         const double &spinParameter, const double &stepLength)
{
    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0.0;
    double waveFunctionPlus = 0.0;
    double waveFunctionCurrent = TrialWaveFunction(r, alpha, beta, omega, spinParameter);

    double stepLengthFraction = 1.0 / stepLength;

    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i, j) += stepLength;
            rMinus(i, j) -= stepLength;
            waveFunctionMinus = TrialWaveFunction(rMinus, alpha, beta, omega, spinParameter);
            waveFunctionPlus = TrialWaveFunction(rPlus, alpha, beta, omega, spinParameter);
            QForce(i, j) = (waveFunctionPlus - waveFunctionMinus);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    QForce = QForce * stepLengthFraction / waveFunctionCurrent;
} */
