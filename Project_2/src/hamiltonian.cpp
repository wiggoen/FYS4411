#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
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
    return 1 / distance;
}


double Hamiltonian::LocalEnergyTwoElectrons(const arma::mat &r, const double &alpha, const double &beta,
                                            const double &omega, const double &spinParameter, const bool UseJastrowFactor,
                                            const bool UseFermionInteraction)
{
    double AlphaOmega = alpha*omega;
    double omegaSquaredHalf = 0.5*omega*omega;

    double r_1Squared = r(0, 0)*r(0, 0) + r(0, 1)*r(0, 1);
    double r_2Squared = r(1, 0)*r(1, 0) + r(1, 1)*r(1, 1);
    double alphaSquared = alpha*alpha;

    double energyWithoutJastrow = 2*AlphaOmega + omegaSquaredHalf*(1 - alphaSquared)*(r_1Squared + r_2Squared);

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
                                const double &omega, const double &spinParameter, bool UseJastrowFactor,
                                const bool UseFermionInteraction)
{
    if (nParticles == 2) {
        return LocalEnergyTwoElectrons(r, alpha, beta, omega, spinParameter, UseJastrowFactor, UseFermionInteraction);
    } else
    {
        double LocalEnergy = LocalEnergyMoreParticles(r,nParticles,beta,omega,spinParameter,UseFermionInteraction);
        if (UseJastrowFactor) {
            Wavefunction wavefunc;
            double slater = wavefunc.SlaterDeterminant(nParticles,r,alpha,omega);
            LocalEnergy *= slater;
            }
        return LocalEnergy;
        }
}


double Hamiltonian::LocalEnergyMoreParticles(const arma::mat &r, const int &nParticles, const double &beta,
                                             const double omega, const double &spinParameter, const bool UseFermionInteraction)
{
    double energy = 0;
    double halfOmegaSquared = 0.5*omega*omega;
    for (int j=0; j<nParticles; j++)
    {
        for (int i=0; i<nParticles; i++)
        {
            double Ri = sqrt(r(i,0)*r(i,0)+r(i,1)*r(i,1));
            double Rj = sqrt(r(j,0)*r(j,0)+r(j,1)*r(j,1));

            double top = -0.5*2*spinParameter*beta;
            double bot = (1+beta*abs(Ri-Rj))*(1+beta*abs(Ri-Rj))*(1+beta*abs(Ri-Rj));
            energy += top/bot + halfOmegaSquared*Ri*Ri;
        }
    }
    if (!UseFermionInteraction)
    /* Without interaction */
    {
        return energy;
    } else {
        double interactionTerm = 0;
        for (int j=0; j<nParticles; j++)
        {
            for (int i=0; i<nParticles; i++)
            {
                double Ri = sqrt(r(i,0)*r(i,0)+r(i,1)*r(i,1));
                double Rj = sqrt(r(j,0)*r(j,0)+r(j,1)*r(j,1));
                interactionTerm += 1.0/abs(Ri-Rj);
            }
        }
        return energy + interactionTerm;
    }
}


double Hamiltonian::NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                         const double &alpha, const double &beta, const double &omega,
                                         const double &spinParameter, const bool UseJastrowFactor,
                                         const bool UseNumericalPotentialEnergy)
{
    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0.0;
    double waveFunctionPlus  = 0.0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, alpha, beta, omega, spinParameter, UseJastrowFactor);

    double h = 1e-4;
    double hSquared = h*h;

    /* Kinetic energy */
    double kineticEnergy = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rPlus(i, j) += h;
            rMinus(i, j) -= h;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            waveFunctionPlus = Wavefunction::TrialWaveFunction(rPlus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            kineticEnergy -= (waveFunctionPlus + waveFunctionMinus - 2.0 * waveFunctionCurrent);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    kineticEnergy = 0.5 * kineticEnergy / (waveFunctionCurrent * hSquared);

    if (!UseNumericalPotentialEnergy)
    {
        return kineticEnergy;
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
        return kineticEnergy + externalPotential;
    }
}
