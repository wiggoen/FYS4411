#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include "inc/hermite.h"
#include <armadillo>


Hamiltonian::Hamiltonian( void )
{

}


Hamiltonian::~Hamiltonian( void )
{

}


double Hamiltonian::ParticleDistance(const arma::rowvec &r_i, const arma::rowvec &r_j)
{
    return arma::norm(r_i - r_j);
}


double Hamiltonian::RepulsiveInteraction(const arma::rowvec &r_i, const arma::rowvec &r_j)
{
    double distance = ParticleDistance(r_i, r_j);
    return 1.0/distance;
}


double Hamiltonian::LocalEnergyTwoElectrons(const arma::mat &r, const double &alpha, const double &beta,
                                            const double &omega, const double &spinParameter, const bool &UseJastrowFactor,
                                            const bool &UseFermionInteraction)
{
    double AlphaOmega = alpha*omega;
    double omegaSquaredHalf = 0.5*omega*omega;

    double r1Squared = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1);
    double r2Squared = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1);
    double alphaSquared = alpha*alpha;

    double energyWithoutJastrow = 2*AlphaOmega + omegaSquaredHalf*(1 - alphaSquared)*(r1Squared + r2Squared);

    if (!UseJastrowFactor)
    {
        /* Energy without Jastrow term */
        if (!UseFermionInteraction) {
            /* Without interaction */
            return energyWithoutJastrow;
        } else
        {
            /* With interaction */
            double interaction = RepulsiveInteraction(r.row(0), r.row(1));
            return energyWithoutJastrow + interaction;
        }
    } else
    {
        /* Energy with Jastrow term */
        double r_12 = arma::norm(r.row(0) - r.row(1));
        double betaR_12 = beta*r_12;
        double denominator = 1 + betaR_12;
        double denominatorSquared = denominator*denominator;

        double parentheses = -AlphaOmega*r_12 + spinParameter/denominatorSquared + (1 - betaR_12)/(r_12*denominator);
        double JastrowTerm = (spinParameter/denominatorSquared) * parentheses;

        if (!UseFermionInteraction) {
            /* Without interaction */
            return energyWithoutJastrow - JastrowTerm;
        } else {
            /* With interaction */
            double interaction = RepulsiveInteraction(r.row(0), r.row(1));;
            return energyWithoutJastrow - JastrowTerm + interaction;
        }
    }
}


double Hamiltonian::LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha, const double &beta,
                                const double &omega, const double &spinParameter, bool &UseJastrowFactor,
                                const bool &UseFermionInteraction)
{
    if (nParticles == 2) {
        return LocalEnergyTwoElectrons(r, alpha, beta, omega, spinParameter, UseJastrowFactor, UseFermionInteraction);
    } else
    {
        double LocalEnergy = LocalEnergyMoreParticles(r, nParticles, beta, omega, spinParameter, UseFermionInteraction);
        if (UseJastrowFactor) {
            arma::mat positions = Wavefunction::Positions(nParticles);
            double slater = SlaterEnergy(r, nParticles, omega, positions);
            LocalEnergy *= slater;
            }
        return LocalEnergy;
        }
}


double Hamiltonian::LocalEnergyMoreParticles(const arma::mat &r, const int &nParticles, const double &beta,
                                             const double &omega, const double &spinParameter, const bool &UseFermionInteraction)
{
    double energy = 0;
    double halfOmegaSquared = 0.5*omega*omega;
    for (int j=0; j<nParticles; j++)
    {
        for (int i=0; i<j; i++)
        {
            double Ri = sqrt(r(i,0)*r(i,0)+r(i,1)*r(i,1));
            double Rj = sqrt(r(j,0)*r(j,0)+r(j,1)*r(j,1));

            double top = -0.5*2*spinParameter*beta;
            double bot = (1+beta*abs(Ri-Rj))*(1+beta*abs(Ri-Rj))*(1+beta*abs(Ri-Rj));
            energy += top/bot + halfOmegaSquared*Ri*Ri;
        }
    }
    double interactionTerm = 0;
    if (!UseFermionInteraction)
    /* Without interaction */
    {
        return energy;
    } else {
        for (int j=0; j<nParticles; j++)
        {
            for (int i=0; i<nParticles; i++)
            {
                double Ri = sqrt(r(i,0)*r(i,0)+r(i,1)*r(i,1));
                double Rj = sqrt(r(j,0)*r(j,0)+r(j,1)*r(j,1));
                interactionTerm += 1.0/abs(Ri-Rj);
            }
        }
    }
    return energy + interactionTerm;
}

double Hamiltonian::SlaterEnergy(const arma::mat &r, const int &nParticles, const double &omega, arma::mat &positions)
{
    double slaterEnergy = 0;
    for (int iParticle = 0; iParticle < nParticles; iParticle++)
    {
        for (int jPosition=0; jPosition < nParticles; jPosition++)
        {
            double nx = positions(jPosition,0);
            double ny = positions(jPosition,1);
            double xPosition = (r(iParticle,0));
            double yPosition = (r(iParticle,1));
            slaterEnergy += DerivativeSlater(omega, xPosition, yPosition, nx, ny);
        }
    }
    //arma::mat D = Wavefunction::SlaterDeterminant(nParticles,r,omega);
    return slaterEnergy;
}

double Hamiltonian::DerivativeSlater(const double &omega, const double &xPosition, const double &yPosition, const int &nx, const int &ny)
{
    double x = xPosition;
    double y = yPosition;
    double normFactor    = 1.0;
    double omegaSquared  = omega*omega;
    double derivativeExp = normFactor*omegaSquared*(0.5+x*x+y*y);
    double firstHermit   = DerivativeHermite(nx,omega,x);
    double secondHermit  = DerivativeHermite(ny,omega,y);
    return firstHermit + secondHermit + derivativeExp;
}

double Hamiltonian::DerivativeHermite(const int &n, const double &omega, const double &x)
{
    if (n==0) {return 0;}
    else if (n==1) {return sqrt(omega);}
    else if (n==2) {return 2*omega*x;}
    else {std::cerr << "Something went wrong in the DerivateHermite function" << std::endl; exit(1);}
}

arma::rowvec Hamiltonian::NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                               const double &alpha, const double &beta, const double &omega,
                                               const double &spinParameter, const bool &UseJastrowFactor,
                                               const bool &UseFermionInteraction, const bool &UseNumericalPotentialEnergy)
{
    arma::rowvec energyVector = arma::zeros<arma::rowvec>(nDimensions+1);

    arma::mat rPlus  = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus  = r;
    rMinus = r;

    double waveFunctionMinus   = 0.0;
    double waveFunctionPlus    = 0.0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, alpha, beta, omega, spinParameter, UseJastrowFactor);

    double h = 1e-4;
    double hSquared = h*h;

    /* Kinetic energy */
    double kineticEnergy = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i, j)  += h;
            rMinus(i, j) -= h;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            waveFunctionPlus  = Wavefunction::TrialWaveFunction(rPlus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            kineticEnergy -= (waveFunctionPlus + waveFunctionMinus - 2.0 * waveFunctionCurrent);
            rPlus(i, j)  = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    kineticEnergy = 0.5 * kineticEnergy / (waveFunctionCurrent * hSquared);

    if (!UseNumericalPotentialEnergy)
    {
        energyVector.at(0) = kineticEnergy;
        energyVector.at(1) = kineticEnergy;
        energyVector.at(2) = 0;
        return energyVector;
    } else
    {
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
            externalPotential += 0.5 * omega*omega * rSquared;
        }
        energyVector.at(0) = kineticEnergy + externalPotential;
        energyVector.at(1) = kineticEnergy;
        energyVector.at(2) = externalPotential;
        if (!UseFermionInteraction)
        {
            /* Without interaction */
            return energyVector;
        } else
        {
            /* With interaction */
            double r_12Squared = (r(0, 0)-r(1, 0))*(r(0, 0)-r(1, 0)) + (r(0, 1)-r(1, 1))*(r(0, 1)-r(1, 1));
            double numericalExternalPotential = 1.0/sqrt(r_12Squared);
            energyVector.at(0) += numericalExternalPotential;
            energyVector.at(2) += numericalExternalPotential;
            return energyVector;
        }
    }
}
