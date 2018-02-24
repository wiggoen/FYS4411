#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include "armadillo"


class Hamiltonian
{
public:
    Hamiltonian(int nParticles, const arma::vec &x);
    ~Hamiltonian();
    double Derivative(int nParticles, double h, arma::vec x);
    double DoubleDerivative(int nParticles, double h, arma::vec x);
    class Wavefunction*  wf = nullptr;

};

#endif // HAMILTONIAN_H
