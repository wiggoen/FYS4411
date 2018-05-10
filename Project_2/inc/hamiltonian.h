#pragma once
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>


class Hamiltonian
{
public:
    Hamiltonian( void );
    ~Hamiltonian( void );
    static double LocalEnergyTwoElectrons(const arma::mat &r, const double &alpha, const double &beta,
                                          const double &omega, const double &spinParameter, bool UseJastrowFactor);
    static double LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha,
                              const double &beta, const double &omega, const double &spinParameter, bool UseJastrowFactor);
    static double NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                       const double &alpha, const double &beta, const double &omega,
                                       const double &spinParameter, const double &stepLength, const bool UseJastrowFactor);
private:
};

#endif /* HAMILTONIAN_H */

