#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>


class Hamiltonian
{
public:
    Hamiltonian( void );
    ~Hamiltonian( void );
    static double LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha, const double &beta, const double &omega, const double &a);
private:
};

#endif /* HAMILTONIAN_H */

