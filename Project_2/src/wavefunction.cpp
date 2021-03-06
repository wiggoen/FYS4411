#include "inc/wavefunction.h"
#include "inc/hamiltonian.h"
#include "inc/hermite.h"


Wavefunction::Wavefunction( void )
{

}


Wavefunction::~Wavefunction( void )
{

}


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


double Wavefunction::JastrowWavefunction(const arma::mat &r, const int &nParticles, const double &beta)
{
    double exponential = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            if (i < j)
            {
                double rij = arma::norm(r.row(i) - r.row(j));
                double spinParameter = Hamiltonian::getSpinParameter(nParticles, i, j);
                exponential += spinParameter*rij/(1 + beta*rij);
            }
        }
    }
    return exp(exponential);
}


double Wavefunction::JastrowRatio(const arma::mat &rNew, const arma::mat &rOld, const int &nParticles, const double &beta)
{
    double jastrowNew = JastrowWavefunction(rNew, nParticles, beta);
    double jastrowOld = JastrowWavefunction(rOld, nParticles, beta);

    return jastrowNew/jastrowOld;
}


double Wavefunction::phi(const arma::mat &r, const double &alpha, const double &omega, const int &nx, const int &ny,
                         const int &k)
{
    /* Single particle states, given by Hermite polynomials */
    double sqrtAlphaOmega = sqrt(alpha*omega);
    double xPosition = r(k, 0);
    double yPosition = r(k, 1);
    double hermiteNx = Hermite::H(nx, sqrtAlphaOmega*xPosition);
    double hermiteNy = Hermite::H(ny, sqrtAlphaOmega*yPosition);
    return hermiteNx*hermiteNy*exp(-0.5*alpha*omega*(xPosition*xPosition + yPosition*yPosition));
}


arma::rowvec Wavefunction::phiGradient(const int &nDimensions, const double &alpha, const double &omega, const double &x,
                                       const double &y, const int &nx, const int &ny)
{
    /* Returns phiGradient */
    arma::rowvec gradient = arma::zeros<arma::rowvec>(nDimensions);

    double alphaOmega = alpha*omega;
    double sqrtAlphaOmega = sqrt(alphaOmega);

    double exponential = exp(-0.5*alphaOmega*(x*x + y*y));

    double H_nx = Hermite::H(nx, sqrtAlphaOmega*x);
    double H_ny = Hermite::H(ny, sqrtAlphaOmega*y);

    double derivativeH_nx = Hermite::DerivativeHermite(nx, sqrtAlphaOmega*x);
    double derivativeH_ny = Hermite::DerivativeHermite(ny, sqrtAlphaOmega*y);

    gradient.at(0) = H_ny*exponential*(derivativeH_nx - H_nx*alphaOmega*x);
    gradient.at(1) = H_nx*exponential*(derivativeH_ny - H_ny*alphaOmega*y);

    return gradient;
}


double Wavefunction::phiLaplace(const double &alpha, const double &omega, const double &x, const double &y,
                                const int &nx, const int &ny)
{
    /* Returns phiLaplace */
    double alphaOmega = alpha*omega;
    double sqrtAlphaOmega = sqrt(alphaOmega);

    double exponential = exp(-0.5*alphaOmega*(x*x + y*y));

    double H_nx = Hermite::H(nx, sqrtAlphaOmega*x);
    double H_ny = Hermite::H(ny, sqrtAlphaOmega*y);

    double derivativeH_nx = Hermite::DerivativeHermite(nx, sqrtAlphaOmega*x);
    double derivativeH_ny = Hermite::DerivativeHermite(ny, sqrtAlphaOmega*y);

    double doubleDerivativeH_nx = Hermite::DoubleDerivativeHermite(nx, sqrtAlphaOmega*x);
    double doubleDerivativeH_ny = Hermite::DoubleDerivativeHermite(ny, sqrtAlphaOmega*y);

    double doubleDerivatives = H_ny*doubleDerivativeH_nx + H_nx*doubleDerivativeH_ny;
    double derivatives = -2.0*alphaOmega*(H_ny*derivativeH_nx*x + H_nx*derivativeH_ny*y);
    double positions = H_nx*H_ny*alphaOmega*(alphaOmega*(x*x + y*y) - 2.0);

    return exponential*(doubleDerivatives + derivatives + positions);
}


void Wavefunction::QuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles, const int &nDimensions,
                                const double &alpha, const double &beta, const double &omega, const bool &UseJastrowFactor,
                                const arma::mat &InverseSlaterUp, const arma::mat &InverseSlaterDown, const int &k)
{
    arma::mat slaterGradient = Hamiltonian::SlaterGradient(r, nParticles, nDimensions, alpha, omega, InverseSlaterUp,
                                                           InverseSlaterDown);

    arma::mat jastrowGradient = arma::zeros<arma::mat>(nParticles, nDimensions);
    if (UseJastrowFactor)
    {
        jastrowGradient = Hamiltonian::JastrowGradient(r, nParticles, nDimensions, beta);
    }

    QForce.row(k) = 2*(slaterGradient.row(k) + jastrowGradient.row(k));
}


double Wavefunction::DerivativePsiOfAlpha(const arma::mat &r, const double &omega)
/* Returns 1/Psi * dPsi/dAlpha */
{
    /* Without Jastrow factor for two electrons */
    double r_1Squared = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1);
    double r_2Squared = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1);
    return -0.5*omega*(r_1Squared + r_2Squared);
}


double Wavefunction::DerivativePsiOfBeta(const arma::mat &r, const double &beta)
/* Returns 1/Psi * dPsi/dBeta */
{
    double spinParameter = 1.0;

    /* With Jastrow factor for two electrons */
    double r_12 = arma::norm(r.row(0) - r.row(1));
    double denominator = (1 + beta*r_12)*(1 + beta*r_12);
    return -(spinParameter*(r_12*r_12))/denominator;
}


double Wavefunction::DerivativePsiManyOfAlpha(const arma::mat &r, const int &nParticles, const double &alpha,
                                              const double &omega)
/* Returns 1/Psi * dPsi/dAlpha */
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


double Wavefunction::DerivativePsiManyOfBeta(const arma::mat &r, const int &nParticles, const double &beta)
/* Returns 1/Psi * dPsi/dBeta */
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
                double spinParameter  = Hamiltonian::getSpinParameter(nParticles, k, j);
                double fraction_kj    = spinParameter * (distanceRkj*distanceRkj) / (denominator_kj*denominator_kj);
                sum -= fraction_kj;
            }
        }
    }
    return sum;
}


double Wavefunction::NumericalTrialWaveFunction(const arma::mat &r, const int &nParticles, const double &alpha,
                                                const double &beta, const double &omega, const bool &UseJastrowFactor)
{
    arma::mat QuantumNumber = Hermite::QuantumNumbers();

    arma::mat slaterMatrixUp   = arma::zeros<arma::mat>(nParticles/2, nParticles/2);
    arma::mat slaterMatrixDown = arma::zeros<arma::mat>(nParticles/2, nParticles/2);

    /* Slater matrix up */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            double phi = Wavefunction::phi(r, alpha, omega, nx, ny, i);
            slaterMatrixUp(j, i) = phi;
        }
    }
    /* Slater matrix down */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            double phi = Wavefunction::phi(r, alpha, omega, nx, ny, i + nParticles/2);
            slaterMatrixDown(j, i) = phi;
        }
    }

    double jastrowFactor = 1.0;
    if (UseJastrowFactor)
    {
        double jastrowExponential = 0.0;
        for (int k = 0; k < nParticles; k++)
        {
            for (int j = 0; j < nParticles; j++)
            {
                if (j != k)
                {
                    double distanceRkj   = Hamiltonian::ParticleDistance(r.row(k), r.row(j));
                    double denominator   = (1 + beta*distanceRkj);
                    double spinParameter = Hamiltonian::getSpinParameter(nParticles, k, j);
                    jastrowExponential  += spinParameter * distanceRkj / (denominator*denominator);
                }
            }
            jastrowFactor = exp(jastrowExponential);
        }

    }

    double slaterDeterminantUp   = arma::det(slaterMatrixUp);
    double slaterDeterminantDown = arma::det(slaterMatrixDown);

    return slaterDeterminantUp * slaterDeterminantDown * jastrowFactor;
}


void Wavefunction::NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                         const int &nDimensions, const double &alpha, const double &beta,
                                         const double &omega, const bool &UseJastrowFactor)
{
    arma::mat rPlus  = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus  = r;
    rMinus = r;

    double waveFunctionMinus   = 0.0;
    double waveFunctionPlus    = 0.0;
    double waveFunctionCurrent = Wavefunction::NumericalTrialWaveFunction(r, nParticles, alpha, beta, omega, UseJastrowFactor);

    double h = 1e-4;

    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i, j)  += h;
            rMinus(i, j) -= h;
            waveFunctionMinus = Wavefunction::NumericalTrialWaveFunction(rMinus, nParticles, alpha, beta, omega, UseJastrowFactor);
            waveFunctionPlus  = Wavefunction::NumericalTrialWaveFunction(rPlus, nParticles, alpha, beta, omega, UseJastrowFactor);
            QForce(i, j) = (waveFunctionPlus - waveFunctionMinus);
            rPlus(i, j)  = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    QForce = QForce / (waveFunctionCurrent * h);
}
