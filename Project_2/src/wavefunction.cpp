#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include "inc/hermite.h"


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
        double r_12 = arma::norm(r.row(0) - r.row(1));
        double jastrow = (spinParameter*r_12)/(1 + beta*r_12);
        return exp(unperturbed + jastrow);
    }
}

double Wavefunction::TrialWaveFunctionManyParticles(const arma::mat &r, const double nParticles, const double &alpha, const double &beta, const double &omega,
                                                    const double &spinParameter, const bool &UseJastrowFactor)
{
    double wavefuncProd = 1;
    double slater       = 1;
    for (int i=0; i<nParticles; i++)
    {
        for (int j=0; j<nParticles; j++)
        {
            double Rij = arma::norm(r.row(i) - r.row(j));
            double expontential = spinParameter*Rij/(1+beta*Rij);
            wavefuncProd *= exp(expontential);
        }
    }
    if (UseJastrowFactor)
    {
        slater = SlaterDeterminant(nParticles, r, alpha, omega);
    }
    return wavefuncProd*slater;
}


double Wavefunction::SlaterDeterminant(const int nParticles, const arma::mat &r, const double &alpha, const double &omega)
{
    /* Create NxN matrix */
    arma::mat slater = arma::zeros<arma::mat>(nParticles, nParticles);

    /* Fill Jastrow matrix */
    int nx=0;
    int ny=0;    //HOW TO FIX THESE??
    for (int iParticle = 0; iParticle < nParticles; iParticle++)
    {
        if (iParticle<2){nx=0;}
        else if (iParticle<6){nx=1;}
        else if (iParticle<12){nx=2;}
        else {std::cout << "Too many particles!" << std::endl
                        << "Please use nParticles=2, 6 or 12." << std::endl;}
        for (int jPosition=0; jPosition < nParticles; jPosition++)
        {
            if (jPosition<2){ny=0;}
            else if (jPosition<6){ny=1;}
            else if (jPosition<12){ny=2;}
            slater(iParticle,jPosition) = phi(iParticle, jPosition, r, nParticles, alpha, omega, nx, ny);
        }
    }
    return det(slater);
}


double Wavefunction::phi(const int &i, const int &j, const arma::mat &r, const int nParticles, const double &alpha, const double &omega,
                         const int &nx, const int &ny)
{
    /* Single particle states, given by Hermite polynomials */
    double sqrtOmega = sqrt(omega);
    double xPosition = r(i,0);
    double yPosition = r(i,1);
    double hermiteNx = H(sqrtOmega*xPosition, nx);
    double hermiteNy = H(sqrtOmega*yPosition, ny);
    return hermiteNx*hermiteNy*exp(-omega*(xPosition*xPosition + yPosition*yPosition)/2.0);
}




void Wavefunction::QuantumForce(const arma::mat &r, arma::mat &QForce, const double &alpha, const double &beta,
                                const double &omega, const double &spinParameter, const bool &UseJastrowFactor)
{
    if (!UseJastrowFactor) {
        /* Without Jastrow factor */
        QForce = -2*alpha*omega*r;
    }
    else {
        /* With Jastrow factor */
        double r_12 = arma::norm(r.row(0) - r.row(1));
        arma::rowvec JastrowDependence = (spinParameter*(r.row(0) - r.row(1)))/(r_12 * ((1 + beta*r_12)*(1 + beta*r_12)));
        QForce.row(0) = 2*(-alpha*omega*r.row(0) + JastrowDependence);
        QForce.row(1) = 2*(-alpha*omega*r.row(1) - JastrowDependence);
    }
}


double Wavefunction::DerivativePsiOfAlpha(const arma::mat &r, const double &omega)
/* Returns 1/Psi * dPsi/dAlpha */
{
    /* Without Jastrow factor */
    double r_1Squared = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1);
    double r_2Squared = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1);
    return -0.5*omega*(r_1Squared + r_2Squared);
}


double Wavefunction::DerivativePsiOfBeta(const arma::mat &r, const double &beta, const double &spinParameter)
/* Returns 1/Psi * dPsi/dBeta */
{
    /* With Jastrow factor */
    double r_12 = arma::norm(r.row(0) - r.row(1));
    double denominator = (1 + beta*r_12)*(1 + beta*r_12);
    return -(spinParameter*(r_12*r_12))/denominator;
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
