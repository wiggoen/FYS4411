#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>


class Hamiltonian
{
public:
    Hamiltonian( void );
    ~Hamiltonian( void );
    static double LocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions, const double &alpha);
    static double NumericalLocalEnergy(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                       const double &alpha, const double &stepLength, const double &beta);
    static double LocalEnergyInteraction(const arma::mat &r, const int &nParticles, const int &nDimensions,
                                         const double &alpha, const double &beta, const double &a);
    static double ParticleDistance(const arma::rowvec &r_k, const arma::rowvec &r_j);
    static arma::rowvec VectorSum(const arma::mat &r, const int &nParticles, const int &nDimensions, const double &a,
                                  const int &k);
private:
    static double DoubleSum(const arma::mat &r, const int &nParticles, const double &a, const int &k);
    static double DerivativeSum(const arma::mat &r, const int &nParticles, const double &a, const int &k);
    static double RepulsivePotential(const arma::mat &r, const int &nParticles, const double &a, const int &k);
};

#endif /* HAMILTONIAN_H */

