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

    static double getSpinParameter(const int &nParticles, const int &i, const int &j);

    static arma::rowvec LocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                    const double &alpha, const double &beta, const double &omega,
                                    const bool &UseJastrowFactor, const bool &UseFermionInteraction,
                                    const arma::mat &InverseSlaterUp, const arma::mat &InverseSlaterDown);

    static arma::mat SlaterGradient(const arma::mat &r, const int &nParticles, const int &nDimensions, const double &alpha,
                                    const double &omega, const arma::mat &InverseSlaterUp,
                                    const arma::mat &InverseSlaterDown);

    static double SlaterLaplacian(const arma::mat &r, const int &nParticles, const double &alpha, const double &omega,
                                  const arma::mat &InverseSlaterUp, const arma::mat &InverseSlaterDown);

    static arma::mat JastrowGradient(const arma::mat &r, const int &nParticles, const int &nDimensions, const double &beta);

    static double JastrowLaplacian(const arma::mat &r, const int &nParticles, const double &beta);

    static arma::rowvec NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                             const double &alpha, const double &beta, const double &omega,
                                             const bool &UseJastrowFactor, const bool &UseFermionInteraction,
                                             const bool &UseNumericalPotentialEnergy);

private:
};

#endif /* HAMILTONIAN_H */

