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


//double Wavefunction::TrialWaveFunctionManyParticles(const arma::mat &r, const int &nParticles, const double &beta,
//                                                    const arma::mat &spinMatrix, const bool &UseJastrowFactor,
//                                                    const arma::mat &SlaterUp, const arma::mat &SlaterDown)
//{
//    double slaterWavefunction = SlaterUp*SlaterDown;
//
//    double jastrowWavefunction = 1.0;
//    if (UseJastrowFactor)
//    {
//        jastrowWavefunction = JastrowWavefunction(r, nParticles, beta, spinMatrix);
//    }
//    return slaterWavefunction*jastrowWavefunction;
//}


double Wavefunction::SlaterRatio(const arma::mat &rNew, const int &nParticles, const double &alpha, const double &omega,
                                 const arma::mat &InverseSlaterUp, const arma::mat &InverseSlaterDown, const int &i)
{
    const arma::mat QuantumNumber = Hermite::QuantumNumbers();
    double slaterRatio = 0.0;

    if (i < nParticles/2)
    {
        /* Electron with spin up moved */
        for (int k = 0; k < nParticles/2; k++)
        {
            int nx = QuantumNumber(k, 0);
            int ny = QuantumNumber(k, 1);
            slaterRatio += phi(rNew, alpha, omega, nx, ny, i)*InverseSlaterUp(k, i);
        }

    } else
    {
        /* Electron with spin down moved */
        for (int k = 0; k < nParticles/2; k++)
        {
            int nx = QuantumNumber(k, 0);
            int ny = QuantumNumber(k, 1);
            slaterRatio += phi(rNew, alpha, omega, nx, ny, i)*InverseSlaterDown(k, i-nParticles/2);
        }
    }
    return slaterRatio;
}


double Wavefunction::JastrowWavefunction(const arma::mat &r, const int &nParticles, const double &beta,
                                         const arma::mat &spinMatrix)
{
    double exponential = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            if (i < j)
            {
                double rij = arma::norm(r.row(i) - r.row(j));
                exponential += spinMatrix(i, j)*rij/(1 + beta*rij);
            }
        }
    }
    return exp(exponential);
}


double Wavefunction::JastrowRatio(const arma::mat &rNew, const arma::mat &rOld, const int &nParticles, const double &beta,
                                  const arma::mat &spinMatrix)
{
    double jastrowNew = JastrowWavefunction(rNew, nParticles, beta, spinMatrix);
    double jastrowOld = JastrowWavefunction(rOld, nParticles, beta, spinMatrix);

    return jastrowNew/jastrowOld;
}


double Wavefunction::phi(const arma::mat &r, const double &alpha, const double &omega, const int &nx, const int &ny, const int &k)
{
    /* Single particle states, given by Hermite polynomials */
    double sqrtAlphaOmega = sqrt(alpha*omega);
    double xPosition = r(k, 0);
    double yPosition = r(k, 1);
    double hermiteNx = Hermite::H(nx, sqrtAlphaOmega*xPosition);
    double hermiteNy = Hermite::H(ny, sqrtAlphaOmega*yPosition);
    return hermiteNx*hermiteNy*exp(-0.5*alpha*omega*(xPosition*xPosition + yPosition*yPosition));
}


arma::rowvec Wavefunction::phiGradient(const int &nDimensions, const double &alpha, const double &omega, const double &x, const double &y,
                                       const int &nx, const int &ny)
{
    /* Returns phiGradient/phi */
    arma::rowvec gradient = arma::zeros<arma::rowvec>(nDimensions);

    double alphaOmega = alpha*omega;
    double sqrtAlphaOmega = sqrt(alphaOmega);

    double inverseH_nx = 1.0/Hermite::H(nx, sqrtAlphaOmega*x);
    double inverseH_ny = 1.0/Hermite::H(ny, sqrtAlphaOmega*y);

    double derivativeH_nx = Hermite::DerivativeHermite(nx, sqrtAlphaOmega*x);
    double derivativeH_ny = Hermite::DerivativeHermite(ny, sqrtAlphaOmega*y);

    gradient.at(0) = inverseH_nx*derivativeH_nx - alphaOmega*x;
    gradient.at(1) = inverseH_ny*derivativeH_ny - alphaOmega*y;

    return gradient;
}


double Wavefunction::phiLaplace(const double &alpha, const double &omega, const double &x, const double &y,
                                const int &nx, const int &ny)
{
    /* Returns phiLaplace/phi */
    double alphaOmega = alpha*omega;
    double sqrtAlphaOmega = sqrt(alphaOmega);

    double inverseH_nx = 1.0/Hermite::H(nx, sqrtAlphaOmega*x);
    double inverseH_ny = 1.0/Hermite::H(ny, sqrtAlphaOmega*y);

    double derivativeH_nx = Hermite::DerivativeHermite(nx, sqrtAlphaOmega*x);
    double derivativeH_ny = Hermite::DerivativeHermite(ny, sqrtAlphaOmega*y);

    double doubleDerivativeH_nx = Hermite::DoubleDerivativeHermite(nx, sqrtAlphaOmega*x);
    double doubleDerivativeH_ny = Hermite::DoubleDerivativeHermite(ny, sqrtAlphaOmega*y);

    double doubleDerivatives = inverseH_nx*doubleDerivativeH_nx + inverseH_ny*doubleDerivativeH_ny;
    double derivatives = -2.0*alphaOmega*(inverseH_nx*derivativeH_nx*x + inverseH_ny*derivativeH_ny*y);
    double positions = alphaOmega*(alphaOmega*(x*x + y*y) - 2.0);

    double laplacian = doubleDerivatives+derivatives+positions;
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


double Wavefunction::DerivativePsiManyOfAlpha(const arma::mat &r, const int &nParticles, const double &alpha, const double &omega)
{
    double sum = 0.0;

    const arma::mat QuantumNumber = Hermite::QuantumNumbers();

    double factor = 0.5*sqrt(omega/alpha);
    double sqrtAlphaOmega = sqrt(alpha*omega);

    /* Electron with spin up moved */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            double x = r(i, 0);
            double y = r(i, 1);
            double inverseH_nx = 1.0/Hermite::H(nx, sqrtAlphaOmega*x);
            double inverseH_ny = 1.0/Hermite::H(ny, sqrtAlphaOmega*y);
            double alphaDerivativeH_nx = Hermite::AlphaDerivativeHermite(nx, x, alpha, omega);
            double alphaDerivativeH_ny = Hermite::AlphaDerivativeHermite(ny, y, alpha, omega);
            sum += factor*x*inverseH_nx*alphaDerivativeH_nx + factor*y*inverseH_ny*alphaDerivativeH_ny;
        }
    }
    /* Electron with spin down moved */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            double x = r(i + nParticles/2, 0);
            double y = r(i + nParticles/2, 1);
            double inverseH_nx = 1.0/Hermite::H(nx, sqrtAlphaOmega*x);
            double inverseH_ny = 1.0/Hermite::H(ny, sqrtAlphaOmega*y);
            double alphaDerivativeH_nx = Hermite::AlphaDerivativeHermite(nx, x, alpha, omega);
            double alphaDerivativeH_ny = Hermite::AlphaDerivativeHermite(ny, y, alpha, omega);
            sum += factor*x*inverseH_nx*alphaDerivativeH_nx + factor*y*inverseH_ny*alphaDerivativeH_ny;
        }
    }
    return  sum - nParticles;
}


double Wavefunction::DerivativePsiManyOfBeta(const arma::mat &r, const int &nParticles, const double &beta,
                                             const arma::mat &spinMatrix)
{
    double sum = 0.0;

    for (int k = 0; k < nParticles; k++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            if (j != k)
            {
                double distanceRkj    = arma::norm(r.row(k) - r.row(j));
                double denominator_kj =  (1 + beta*distanceRkj);
                double fraction_kj    = spinMatrix(k, j) * (distanceRkj*distanceRkj) / (denominator_kj*denominator_kj);
                sum -= fraction_kj;
            }
        }
    }
    return sum;
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
