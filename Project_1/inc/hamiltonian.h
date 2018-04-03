#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>


class Hamiltonian
{
public:
    Hamiltonian();
    ~Hamiltonian();
    static double LocalEnergy(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha);
    static double NumericalLocalEnergy(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha, double &stepLength);
};

#endif // HAMILTONIAN_H

