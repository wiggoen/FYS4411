#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include "inc/hermite.h"
#include <armadillo>
#include <iostream>

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


double Hamiltonian::getSpinParameter(const int &nParticles, const int &i, const int &j)
{
    if       (i < nParticles/2  && j < nParticles/2)  { return 1.0/3.0; }
    else if  (i >= nParticles/2 && j >= nParticles/2) { return 1.0/3.0; }
    else                                              { return 1.0; }
}


double Hamiltonian::LocalEnergyTwoParticles(const arma::mat &r, const double &alpha, const double &beta,
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


arma::rowvec Hamiltonian::LocalEnergyMoreParticles(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                                   const double &alpha, const double &beta, const double &omega,
                                                   arma::mat &spinMatrix, const bool &UseJastrowFactor,
                                                   const bool &UseFermionInteraction, const arma::mat &InverseSlaterUp,
                                                   const arma::mat &InverseSlaterDown)
{
    arma::rowvec energyVector = arma::zeros<arma::rowvec>(nDimensions+1);

    /* Kinetic energy */
    double slaterLaplacian = SlaterLaplacian(r, nParticles, alpha, omega, InverseSlaterUp, InverseSlaterDown);

    double jastrowLaplacian = 0.0;
    double gradients = 0.0;
    if (UseJastrowFactor)
    {
        jastrowLaplacian = JastrowLaplacian(r, nParticles, beta, spinMatrix);
        arma::mat slaterGradient = SlaterGradient(r, nParticles, nDimensions, alpha, omega, InverseSlaterUp, InverseSlaterDown);
        arma::mat jastrowGradient = JastrowGradient(r, nParticles, nDimensions, beta, spinMatrix);
        for (int i = 0; i < nParticles; i++)
        {
            gradients += arma::dot(slaterGradient.row(i), jastrowGradient.row(i));
        }
    }
    double kineticEnergy = - 0.5*slaterLaplacian - 0.5*jastrowLaplacian - gradients;
    //std::cout << "Ks " << -0.5*slaterLaplacian << "   Kj " << -0.5*jastrowLaplacian << "   Kg " <<  - gradients << std::endl;

    /* External potential */
    double rSquared = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            rSquared += r(i, j) * r(i, j);
        }
    }
    double externalPotential = 0.5 * omega*omega * rSquared;


    /* Repulsive potential (interaction) */
    double repulsivePotential = 0.0;
    if (UseFermionInteraction)
    {
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = i + 1; j < nParticles; j++)
            {
                repulsivePotential += RepulsiveInteraction(r.row(i), r.row(j));
            }
        }
    }
    energyVector.at(0) = kineticEnergy + externalPotential + repulsivePotential;
    energyVector.at(1) = kineticEnergy;
    energyVector.at(2) = externalPotential + repulsivePotential;
    //std::cout << "K  " << kineticEnergy << "   Ve " << externalPotential << "   Vi " << repulsivePotential << std::endl;
    //std::cout << std::endl;
    return energyVector;
}


arma::mat Hamiltonian::SlaterGradient(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                      const double &alpha, const double &omega,
                                      const arma::mat &InverseSlaterUp, const arma::mat &InverseSlaterDown)
{
    arma::mat slaterGradient = arma::zeros<arma::mat>(nParticles, nDimensions);

    arma::mat QuantumNumber = Hermite::QuantumNumbers();

    /* Gradient up */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            double x = r(i, 0);
            double y = r(i, 1);
            arma::rowvec phiGradient = Wavefunction::phiGradient(nDimensions, alpha, omega, x, y, nx, ny);
            slaterGradient.row(i) += phiGradient*InverseSlaterUp(j, i);
        }
    }
    /* Gradient down */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            double x = r(i + nParticles/2, 0);
            double y = r(i + nParticles/2, 1);
            arma::rowvec phiGradient = Wavefunction::phiGradient(nDimensions, alpha, omega, x, y, nx, ny);
            slaterGradient.row(i+nParticles/2) += phiGradient*InverseSlaterDown(j, i);
        }
    }
    return slaterGradient;
}


double Hamiltonian::SlaterLaplacian(const arma::mat &r, const int &nParticles, const double &alpha, const double &omega,
                                    const arma::mat &InverseSlaterUp, const arma::mat &InverseSlaterDown)
{
    arma::mat QuantumNumber = Hermite::QuantumNumbers();
    double laplacianUp   = 0.0;
    double laplacianDown = 0.0;

    /* Laplacian up */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            double x = r(i, 0);
            double y = r(i, 1);
            double phiLaplace = Wavefunction::phiLaplace(alpha, omega, x, y, nx, ny);
            laplacianUp += phiLaplace*InverseSlaterUp(j, i);
        }
    }
    /* Laplacian down */
    for (int i = 0; i < nParticles/2; i++)
    {
        for (int j = 0; j < nParticles/2; j++)
        {
            int nx = QuantumNumber(j, 0);
            int ny = QuantumNumber(j, 1);
            double x = r(i + nParticles/2, 0);
            double y = r(i + nParticles/2, 1);
            double phiLaplace = Wavefunction::phiLaplace(alpha, omega, x, y, nx, ny);
            laplacianDown += phiLaplace*InverseSlaterDown(j, i);
        }
    }
    return laplacianUp + laplacianDown;
}


arma::mat Hamiltonian::JastrowGradient(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                       const double &beta, const arma::mat &spinMatrix)
{
    arma::mat jastrowGradient = arma::zeros<arma::mat>(nParticles, nDimensions);

    for (int k = 0; k < nParticles; k++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            if (j != k)
            {
                double distanceRkj = ParticleDistance(r.row(k), r.row(j));
                double denominator = (1 + beta*distanceRkj);
                double spinParameter = getSpinParameter(nParticles, k, j);
                double fraction = spinParameter / (denominator*denominator);
                jastrowGradient.row(k) += ( ((r.row(k) - r.row(j)) / distanceRkj ) * fraction );
            }
        }
    }
    return jastrowGradient;
}


double Hamiltonian::JastrowLaplacian(const arma::mat &r, const int &nParticles, const double &beta,
                                     const arma::mat &spinMatrix)
{
    double singleSum     = 0.0;
    double doubleSum     = 0.0;
    double spinParameter = 0.0;

    for (int k = 0; k < nParticles; k++)
    {
        for (int j = 0; j < nParticles; j++)
        {
            if (j != k)
            {
                arma::rowvec Rkj      = (r.row(k) - r.row(j));
                double distanceRkj    = arma::norm(Rkj);
                double denominator_kj =  (1 + beta*distanceRkj);
                spinParameter         = getSpinParameter(nParticles, k, j);
                double fraction_kj    = spinParameter / (denominator_kj*denominator_kj);
                singleSum += (fraction_kj/distanceRkj - 2*beta*fraction_kj / denominator_kj);

                for (int i = 0; i < nParticles; i++)
                {
                    if (i != k)
                    {
                        arma::rowvec Rki      = (r.row(k) - r.row(i));
                        double distanceRki    = arma::norm(Rki);
                        double denominator_ki =  (1 + beta*distanceRki);
                        spinParameter         = getSpinParameter(nParticles, k, j);
                        double fraction_ki    = spinParameter / (denominator_ki*denominator_ki);
                        doubleSum += ( (arma::dot(Rki, Rkj) / (distanceRki*distanceRkj) ) * fraction_ki * fraction_kj );
                    }
                }
            }
        }
    }
    return doubleSum + singleSum;
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
