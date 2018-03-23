#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>


class Wavefunction
{
public:
    Wavefunction();
    ~Wavefunction();
    static double TrialWaveFunction(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha);
    static void QuantumForce(const arma::mat &r, arma::mat &QForce, double &alpha);
    static double derivativePsi(const arma::mat &R, int nParticles, int nDimensions, double beta);

};

#endif // WAVEFUNCTION_H
