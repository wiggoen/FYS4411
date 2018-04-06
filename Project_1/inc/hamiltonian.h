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
    static double LocalEnergyInteraction(const arma::mat &r, int &nParticles, int &nDimensions, double &alpha, double &beta);
    static double ParticleDistance(const arma::rowvec &r_k, const arma::rowvec &r_j);
    static arma::rowvec VectorSum(const arma::mat &r, int &nParticles, int &nDimensions, const double a, const int &k);
    static double DerivativeSum(const arma::mat &r, int &nParticles, const double a, const int &k);
    static double Correlation(const arma::mat &r, int &nParticles, double &a, int &k);
    static double RepulsivePotential(const arma::mat &r, int &nParticles, double &a, int &k);
};

#endif // HAMILTONIAN_H

