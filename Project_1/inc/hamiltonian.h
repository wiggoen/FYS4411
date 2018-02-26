#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include "armadillo"


class Hamiltonian
{
public:
    Hamiltonian();
    ~Hamiltonian();
    double localEnergy(const arma::mat &r);

};

#endif // HAMILTONIAN_H

