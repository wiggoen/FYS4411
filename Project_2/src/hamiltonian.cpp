#include "inc/hamiltonian.h"
#include "inc/wavefunction.h"
#include <armadillo>


Hamiltonian::Hamiltonian( void )
{

}


Hamiltonian::~Hamiltonian( void )
{

}


double Hamiltonian::LocalEnergyTwoElectrons(const arma::mat &r, const double &alpha, const double &beta,
                                            const double &omega, const double &spinParameter, bool UseJastrowFactor)
{
    double AlphaOmega = alpha*omega;
    double omegaSquaredHalf = 0.5*omega*omega;

    double r_1Squared = r(0,0)*r(0,0) + r(0,1)*r(0,1);
    double r_2Squared = r(1,0)*r(1,0) + r(1,1)*r(1,1);

    if (!UseJastrowFactor) {
        return 2*AlphaOmega + omegaSquaredHalf*(r_1Squared + r_2Squared);
    } else {
        double r_12 = arma::norm(r.row(0) - r.row(1));
        double betaR_12 = beta*r_12;
        double alphaSquared = alpha*alpha;
        double denominator = (1 + beta*r_12);
        double denominatorSquared = denominator*denominator;

        double parentheses = -AlphaOmega*r_12 + 2/r_12 + spinParameter/denominatorSquared - ((1 + 3*betaR_12))/(r_12*denominator);
        double fractions = spinParameter/denominatorSquared * parentheses;

        return 2*AlphaOmega + (1 - alphaSquared)*omegaSquaredHalf*(r_1Squared + r_2Squared) - fractions;
    }
}


double Hamiltonian::LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha, const double &beta,
                                const double &omega, const double &spinParameter, bool UseJastrowFactor)
{
    if (nParticles == 2) {
        return LocalEnergyTwoElectrons(r, alpha, beta, omega, spinParameter, UseJastrowFactor);
    }
    /* Husk å legg til med / uten Jastrow for flere partikler */
    //if (!UseJastrowFactor) {}
    else {
        double localEnergy = 0;
        double AlphaOmega = alpha*omega;
        double omegaSquaredHalf = 0.5*omega*omega;

        double betaR_ij =  beta*arma::norm(r.row(0) - r.row(1));
        //std::cout << "norm = " << betaR_ij << std::endl;

        for (int i = 0; i < nParticles; i++)
        {
            double numerator = spinParameter*(1 + betaR_ij) - spinParameter*betaR_ij;
            double denominator = (1 + betaR_ij)*(1 + betaR_ij);
            double firstFraction = numerator/denominator;
            //double firstTerm = -0.5*(-AlphaOmega*arma::norm(r.row(i)) + firstFraction);
            double firstTerm = -AlphaOmega*arma::norm(r.row(i)) + firstFraction;

            //double secondTerm = 1 - AlphaOmega*arma::dot(r.row(i), r.row(i)) + firstFraction*arma::norm(r.row(i));

            double secondTerm = firstFraction * (2*beta*arma::norm(r.row(i)) + 2*beta*beta*arma::norm(r.row(i)))/denominator - AlphaOmega;

            //double thirdTerm = arma::norm(r.row(i))*(-AlphaOmega + (((1 + beta - a*beta)*(1 + betaR_ij) - numerator)*2*beta)/(denominator*(1 + betaR_ij)));

            double thirdTerm = omegaSquaredHalf*arma::dot(r.row(i), r.row(i));

            //double fourthTerm = omegaSquaredHalf*arma::dot(r.row(i), r.row(i));

            localEnergy += -0.5*firstTerm*firstTerm * (-0.5)*secondTerm * thirdTerm;
        }
        return localEnergy;
    }
}

double Hamiltonian::NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                         const double &alpha, const double &beta, const double &omega,
                                         const double &spinParameter, const double &stepLength, const bool UseJastrowFactor)
{
    arma::mat rPlus = arma::zeros<arma::mat>(nParticles, nDimensions);
    arma::mat rMinus = arma::zeros<arma::mat>(nParticles, nDimensions);

    rPlus = r;
    rMinus = r;

    double waveFunctionMinus = 0.0;
    double waveFunctionPlus  = 0.0;
    double waveFunctionCurrent = Wavefunction::TrialWaveFunction(r, alpha, beta, omega, spinParameter, UseJastrowFactor);

    //std::cout << waveFunctionCurrent << std::endl;

    double stepLengthSquaredFraction = 1.0 / (stepLength * stepLength);

    /* Kinetic energy */
    double kineticEnergy = 0.0;
    for (int i = 0; i < nParticles; i++) {
        for (int j = 0; j < nDimensions; j++) {
            rPlus(i, j) += stepLength;
            rMinus(i, j) -= stepLength;
            waveFunctionMinus = Wavefunction::TrialWaveFunction(rMinus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            waveFunctionPlus = Wavefunction::TrialWaveFunction(rPlus, alpha, beta, omega, spinParameter, UseJastrowFactor);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2.0 * waveFunctionCurrent);
            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }
    kineticEnergy = 0.5 * stepLengthSquaredFraction * kineticEnergy / waveFunctionCurrent;

    /* External potential */
    /*
    double externalPotential = 0.0;
    double rSquared = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        rSquared = 0.0;
        for (int j = 0; j < nDimensions; j++)
        {
            rSquared += r(i, j) * r(i, j);
        }
        externalPotential += 0.5 * rSquared;
    }*/
    return kineticEnergy;// + externalPotential;
}
