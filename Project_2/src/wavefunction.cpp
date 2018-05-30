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
    /* Two electrons */
    double r_1Squared  = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1);
    double r_2Squared  = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1);
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


double Wavefunction::TrialWaveFunctionManyParticles(const arma::mat &r, const double &nParticles, const double &beta,
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
        //slater = SlaterDeterminant(r, nParticles, omega);
    }
    return wavefuncProd*slater;
}


arma::mat Wavefunction::Positions(const int &nParticles)
{
    arma::mat positions = arma::mat(nParticles/2, 2);
    arma::mat possibleQuantumNumbers = arma::mat(nParticles/3, 1);
    for (int i = 0; i < nParticles/3; i++)                              // Hvorfor /3 ?
    {
        possibleQuantumNumbers(i) = i;
    }
    int nx = 0; int ny = 0;
    for (int i = 0; i < nParticles/2; i++)
    {
        if (i == 1) { nx += 1; } if (i == 2){ nx -= 1; ny += 1; }
        if (i == 3) { nx += 1; } if (i == 4){ nx += 1; ny -= 1; } // CAN THIS BE OPTIMIZED?
        if (i == 5) { nx -= 2; ny += 2; }
        positions(i, 0) = possibleQuantumNumbers(nx);
        positions(i, 1) = possibleQuantumNumbers(ny);
    }
    //std::cout << "Position vector: " << std::endl;
    //std::cout << positions << std::endl << std::endl;
    return positions;
}


int factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


arma::mat Wavefunction::SlaterDeterminant(const arma::mat &r, const int &nParticles, const double &alpha, const double &omega)
{
    /* Create NxN matrix */
    arma::mat positions = Positions(nParticles);
    arma::mat slater = arma::zeros<arma::mat>(nParticles/2, nParticles/2);
    int nx=0; int ny=0;
    double normFactor = 1.0/sqrt(factorial(nParticles/2));
    /*   Fill Jastrow matrix   */
    /* Particle i at position j */
    for (int iParticle = 0; iParticle < nParticles/2; iParticle++)
    {
        for (int jPosition=0; jPosition < nParticles/2; jPosition++)
        {
            nx = positions(jPosition,0);
            ny = positions(jPosition,1);
            slater(iParticle,jPosition) = phi(r, alpha, omega, nx, ny, iParticle);
        }
    }
    return arma::inv(slater)*normFactor;
}


double Wavefunction::phi(const arma::mat &r, const double alpha, const double &omega, const int &nx, const int &ny, const int &j)
{
    /* Single particle states, given by Hermite polynomials */
    double sqrtAlphaOmega = sqrt(alpha*omega);
    double xPosition = r(j, 0);
    double yPosition = r(j, 1);
    double hermiteNx = H(sqrtAlphaOmega*xPosition, nx);
    double hermiteNy = H(sqrtAlphaOmega*yPosition, ny);
    //std::cout << hermiteNx << std::endl;
    return hermiteNx*hermiteNy*exp(-0.5*omega*(xPosition*xPosition + yPosition*yPosition));
}


arma::mat Wavefunction::phiGradient(const double &nParticles, const double &alpha, const double &omega, const double &x, const double &y,
                                    const int &nx, const double &ny)
{
    arma::mat gradient = arma::mat(nParticles, 2);

    double alphaOmega = alpha*omega;
    double sqrtAlphaOmega = sqrt(alphaOmega);

    double exponential = exp(-0.5*alphaOmega*(x*x + y*y));

    double H_nx = H(sqrtAlphaOmega*x, nx);
    double H_ny = H(sqrtAlphaOmega*y, ny);

    double derivativeH_nx = Hamiltonian::DerivativeHermite(nx, sqrtAlphaOmega*x);
    double derivativeH_ny = Hamiltonian::DerivativeHermite(ny, sqrtAlphaOmega*y);

    double derivativePsi_x = H_ny*exponential*(derivativeH_nx - H_nx*alphaOmega*x);
    double derivativePsi_y = H_nx*exponential*(derivativeH_ny - H_ny*alphaOmega*y);

    gradient(0, 0) = derivativePsi_x;
    gradient(0, 1) = derivativePsi_y;
    return gradient;
}


double Wavefunction::phiLaplace(const double &alpha, const double &omega, const double &x, const double &y,
                                const int &nx, const double &ny)
{
    double alphaOmega = alpha*omega;
    double sqrtAlphaOmega = sqrt(alphaOmega);
    double exponential = exp(-0.5*alphaOmega*(x*x + y*y));

    double H_nx = H(sqrtAlphaOmega*x, nx);
    double H_ny = H(sqrtAlphaOmega*y, ny);

    double derivativeH_nx = Hamiltonian::DerivativeHermite(nx, sqrtAlphaOmega*x);
    double derivativeH_ny = Hamiltonian::DerivativeHermite(ny, sqrtAlphaOmega*y);

    double doubleDerivativeH_nx = Hamiltonian::DoubleDerivativeHermite(nx, sqrtAlphaOmega*x);
    double doubleDerivativeH_ny = Hamiltonian::DoubleDerivativeHermite(ny, sqrtAlphaOmega*y);

    double term1 = H_ny*doubleDerivativeH_nx + H_nx*doubleDerivativeH_ny;
    double term2 = -2.0*alphaOmega*(H_ny*derivativeH_nx*x + H_nx*derivativeH_ny*y);
    double term3 = H_nx*H_ny*alphaOmega*(alphaOmega*(x*x + y*y) - 2.0);

    double laplacian = exponential*(term1 + term2 + term3);
    return laplacian;
}


void Wavefunction::QuantumForce(const arma::mat &r, arma::mat &QForce, const double &alpha, const double &beta,
                                const double &omega, const double &spinParameter, const bool &UseJastrowFactor)
{
    if (!UseJastrowFactor) {
        /* Without Jastrow factor */
        QForce = -2*alpha*omega*r;
    }
    else {
        /* With Jastrow factor for two electrons */
        double r_12 = arma::norm(r.row(0) - r.row(1));
        arma::rowvec JastrowDependence = (spinParameter*(r.row(0) - r.row(1)))/(r_12 * ((1 + beta*r_12)*(1 + beta*r_12)));
        QForce.row(0) = 2*(-alpha*omega*r.row(0) + JastrowDependence);
        QForce.row(1) = 2*(-alpha*omega*r.row(1) - JastrowDependence);
    }
}


double Wavefunction::DerivativePsiOfAlpha(const arma::mat &r, const double &omega)
/* Returns 1/Psi * dPsi/dAlpha */
{
    /* Without Jastrow factor for two electrons */
    double r_1Squared = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1);
    double r_2Squared = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1);
    return -0.5*omega*(r_1Squared + r_2Squared);
}


double Wavefunction::DerivativePsiOfBeta(const arma::mat &r, const double &beta, const double &spinParameter)
/* Returns 1/Psi * dPsi/dBeta */
{
    /* With Jastrow factor for two electrons */
    double r_12 = arma::norm(r.row(0) - r.row(1));
    double denominator = (1 + beta*r_12)*(1 + beta*r_12);
    return -(spinParameter*(r_12*r_12))/denominator;
}


void Wavefunction::NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                         const int &nDimensions, const double &alpha, const double &beta,
                                         const double &omega, const double &spinParameter, const bool &UseJastrowFactor)
{
    arma::mat rPlus  = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus  = r;
    rMinus = r;

    double waveFunctionMinus   = 0.0;
    double waveFunctionPlus    = 0.0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, alpha, beta, omega, spinParameter, UseJastrowFactor);

    double h = 1e-4;

    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i, j)  += h;
            rMinus(i, j) -= h;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            waveFunctionPlus  = Wavefunction::TrialWaveFunction(rPlus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            QForce(i, j) = (waveFunctionPlus - waveFunctionMinus);
            rPlus(i, j)  = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    QForce = QForce / (waveFunctionCurrent * h);
}
