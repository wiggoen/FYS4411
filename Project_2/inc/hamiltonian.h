#pragma once
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>


class Hamiltonian
{
public:
    Hamiltonian( void );
    ~Hamiltonian( void );
    static double ParticleDistance(const arma::rowvec &r_i, const arma::rowvec &r_j);
    static double RepulsiveInteraction(const arma::rowvec &r_i, const arma::rowvec &r_j);
    static double LocalEnergyTwoElectrons(const arma::mat &r, const double &alpha, const double &beta,
                                          const double &omega, const double &spinParameter, const bool UseJastrowFactor,
                                          const bool UseFermionInteraction);
    static double LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha,
                              const double &beta, const double &omega, const double &spinParameter, bool UseJastrowFactor,
                              const bool UseFermionInteraction);
    static double NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                       const double &alpha, const double &beta, const double &omega,
                                       const double &spinParameter, const bool UseJastrowFactor,
                                       const bool UseNumericalPotentialEnergy);
    static double LocalEnergyMoreParticles(const arma::mat &r, const int &nParticles, const double &beta,
                                           const double &omega, const double &spinParameter, const bool UseFermionInteraction);
private:
};

#endif /* HAMILTONIAN_H */

