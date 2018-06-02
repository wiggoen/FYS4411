#pragma once
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>


class Wavefunction
{
public:
    Wavefunction( void );
    ~Wavefunction( void );
    static double TrialWaveFunction(const arma::mat &r, const double &alpha, const double &beta,
                                    const double &omega, const double &spinParameter, const bool &UseJastrowFactor);
    static void QuantumForce(const arma::mat &r, arma::mat &QForce, const double &alpha, const double &beta,
                             const double &omega, const double &spinParameter, const bool &UseJastrowFactor);
    static double DerivativePsiOfAlpha(const arma::mat &r, const double &omega);
    static double DerivativePsiOfBeta(const arma::mat &r, const double &beta, const double &spinParameter);
    static double DerivativePsiManyOfAlpha(const double &alpha, const double &omega, const int &nx, const int &ny,
                                                  const double &x, const double &y);
    static double DerivativePsiManyOfBeta(const arma::mat &r, const int &nParticles, const double &beta, const double &omega, const arma::mat &spinMatrix);
    static void NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                      const int &nDimensions, const double &alpha, const double &beta,
                                      const double &omega, const double &spinParameter, const bool &UseJastrowFactor);
    static double TrialWaveFunctionManyParticles(const arma::mat &r, const int &nParticles,
                                                 const double &beta, const arma::mat &spinMatrix,
                                                 const bool &UseJastrowFactor, const arma::mat &SlaterUp,
                                                 const arma::mat &SlaterDown);
    static double phiLaplace(const double &alpha, const double &omega, const double &x, const double &y, const int &nx, const int &ny);
    static arma::rowvec phiGradient(const int &nDimensions, const double &alpha, const double &omega, const double &x, const double &y, const int &nx, const int &ny);

    static double JastrowWavefunction(const arma::mat &r, const int &nParticles, const double &beta, const arma::mat &spinMatrix);

    static double JastrowRatio(const arma::mat &rNew, const arma::mat &rOld, const int &nParticles, const double &beta,
                               const arma::mat &spinMatrix);


    static double phi(const arma::mat &r, const double &alpha, const double &omega, const int &nx, const int &ny, const int &k);

    static double SlaterRatio(const arma::mat &rNew, const int &nParticles, const double &alpha, const double &omega,
                              const arma::mat &InverseSlaterUp, const arma::mat &InverseSlaterDown, const int &i);
};


#endif /* WAVEFUNCTION_H */
