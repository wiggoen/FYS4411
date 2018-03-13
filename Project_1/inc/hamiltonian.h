#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>


class Hamiltonian
{
public:
    Hamiltonian();
    ~Hamiltonian();
    static double LocalEnergy(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha);
};

#endif // HAMILTONIAN_H

