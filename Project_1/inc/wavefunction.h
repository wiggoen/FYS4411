#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>


class Wavefunction
{
public:
    Wavefunction();
    ~Wavefunction();
    static double TrialWaveFunction(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha, double &beta);
    static double TrialWaveFunctionInteraction(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha, double &beta, double &a);
    static void QuantumForce(const arma::mat &r, arma::mat &QForce, double &alpha);
    static void NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, int &nParticles, int &nDimensions, double &alpha, double &stepLength);
    static void QuantumForceInteraction(const arma::mat &r, arma::mat &QForce, double &alpha, int &nParticles, int &nDimensions, double &a, int &i);
    static double DerivativePsi(const arma::mat &r, int &nParticles, int &nDimensions, double &beta);

};

#endif // WAVEFUNCTION_H
