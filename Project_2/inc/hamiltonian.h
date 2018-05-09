#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>


class Hamiltonian
{
public:
    Hamiltonian( void );
    ~Hamiltonian( void );
    static double LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha, const double &beta, const double &omega, const double &a);
    static double NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                             const double &alpha, const double &stepLength, const double &beta, const double &omega, const double &a, const double &constant);
private:
};

#endif /* HAMILTONIAN_H */

