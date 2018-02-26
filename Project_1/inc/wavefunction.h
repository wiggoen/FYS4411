#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>


class Wavefunction
{
public:
    Wavefunction();
    ~Wavefunction();
    double trialWaveFunction(const arma::mat &r, int nParticles, int nDimensions, double alpha);
    void QuantumForce(const arma::mat &r, arma::mat &QForce, double alpha);
};

#endif // WAVEFUNCTION_H
