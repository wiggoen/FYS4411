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
    static void NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                      const int &nDimensions, const double &alpha, const double &beta,
                                      const double &omega, const double &spinParameter, const bool UseJastrowFactor);
    static double TrialWaveFunctionManyParticles(const arma::mat &r, const double nParticles, const double &alpha, const double &beta, const double &omega, const double &spinParameter, const bool &UseJastrowFactor);
    static double SlaterDeterminant(const int nParticles, const arma::mat &r, const double &alpha, const double &omega);
    static double findPossibleNxNy(unsigned int i, const arma::vec &p);
    static double phi(const int &i, const int &j, const arma::mat &r, const int nParticles, const double &alpha, const double &omega, const int &nx, const int &ny);
};


#endif /* WAVEFUNCTION_H */
