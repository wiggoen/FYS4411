#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>


class Wavefunction
{
public:
    Wavefunction( void );
    ~Wavefunction( void );
    static double TrialWaveFunction(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                    const double &alpha, const double &beta);
    static double TrialWaveFunctionInteraction(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                               const double &alpha, const double &beta, const double &a);
    static void QuantumForce(const arma::mat &r, arma::mat &QForce, const double &alpha);
    static void NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                      const int &nDimensions, const double &alpha, const double &stepLength,
                                      const double &beta);
    static void QuantumForceInteraction(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                        const int &nDimensions, const double &alpha, const double &beta, const double &a,
                                        const int &k);
    static double DerivativePsi(const arma::mat &r, const int &nParticles, const int &nDimensions, const double &beta);
};

#endif /* WAVEFUNCTION_H */
