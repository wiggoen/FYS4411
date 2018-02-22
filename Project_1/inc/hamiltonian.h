#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include "armadillo"


class Hamiltonian
{
public:
    Hamiltonian(int nParticles, const arma::mat &);
    ~Hamiltonian();
};

#endif // HAMILTONIAN_H
