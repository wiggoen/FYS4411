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
                                          const double &omega, const double &spinParameter, const bool &UseJastrowFactor,
                                          const bool &UseFermionInteraction);
    static double LocalEnergy(const arma::mat &r, const int &nParticles, const double &alpha,
                              const double &beta, const double &omega, const double &spinParameter, bool &UseJastrowFactor,
                              const bool &UseFermionInteraction);
    static arma::rowvec NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                             const double &alpha, const double &beta, const double &omega,
                                             const double &spinParameter, const bool &UseJastrowFactor,
                                             const bool &UseFermionInteraction, const bool &UseNumericalPotentialEnergy);
    static double LocalEnergyMoreParticles(const arma::mat &r, const int &nParticles, const double &alpha, const double &beta,
                                           const double &omega, const double &spinParameter, const bool &UseFermionInteraction,
                                           const arma::mat InverseSlaterUp, const arma::mat InverseSlaterDown, const int &k);
    static arma::mat GradientSlater(const arma::mat &r, const double &nParticles, const double &alpha, const double &omega,
                                    const double &xPosition, const double &yPosition, const int &nx, const int &ny, const int &iParticle);
    static double LaplaceSlater(const arma::mat &r, const double &nParticles, const double &alpha, const double &omega,
                                const double &xPosition, const double &yPosition, const int &nx, const int &ny, const int &i, const int &k);
    static double JastrowMoreParticles(arma::mat &r, const double &nParticles, const double &spinParameter, const double &beta);
    static double Laplacian(const arma::mat &r, const int &nParticles, const double &alpha, const double &omega,
                            const arma::mat InverseSlaterUp, const arma::mat InverseSlaterDown, const int &k);

    //static double SlaterEnergy(const arma::mat &r, const int &nParticles, const double &alpha, const double &omega, arma::mat &positions);
private:
};

#endif /* HAMILTONIAN_H */

