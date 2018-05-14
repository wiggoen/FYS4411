#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"


Wavefunction::Wavefunction( void )
{

}


Wavefunction::~Wavefunction( void )
{

}


double Wavefunction::TrialWaveFunction(const arma::mat &r, const double &alpha, const double &beta, const double &omega,
                                       const double &spinParameter, const bool &UseJastrowFactor)
{
    double r_1Squared = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1);
    double r_2Squared = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1);
    double unperturbed = -0.5*alpha*omega*(r_1Squared + r_2Squared);

    if (!UseJastrowFactor) {
        /* Without Jastrow factor */
        return exp(unperturbed);
    } else {
        /* With Jastrow factor */
        double r_12 = arma::norm(r.row(0) - r.row(1));;
        double jastrow = (spinParameter*r_12)/(1 + beta*r_12);
        return exp(unperturbed + jastrow);
    }
}


void Wavefunction::QuantumForce(const arma::mat &r, arma::mat &QForce, const double &alpha, const double &omega)
{
    /* Without Jastrow factor */
    //if (!UseJastrowFactor) {

    QForce = -2*alpha*omega*r;

    //}
    //else {
    /* With Jastrow factor */
    //QForce = ?;
    //}
}


double Wavefunction::DerivativePsi(const arma::mat &r, const double &alpha, const double omega)
/* Returns 1/psi * psi' */
{
    /* Without Jastrow */
    /*
    arma::rowvec r_1 = r.row(0);
    arma::rowvec r_2 = r.row(1);

    return -alpha*omega*r;
    */
    return 0;
}


void Wavefunction::NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                         const int &nDimensions, const double &alpha, const double &beta,
                                         const double &omega, const double &spinParameter, const bool UseJastrowFactor)
{
    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0.0;
    double waveFunctionPlus  = 0.0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, alpha, beta, omega, spinParameter, UseJastrowFactor);

    double h = 1e-4;

    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i, j) += h;
            rMinus(i, j) -= h;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            waveFunctionPlus = Wavefunction::TrialWaveFunction(rPlus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            QForce(i, j) = (waveFunctionPlus - waveFunctionMinus);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    QForce = QForce / (waveFunctionCurrent * h);
}
