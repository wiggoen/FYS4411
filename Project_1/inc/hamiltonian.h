#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <armadillo>


class Hamiltonian
{
public:
    Hamiltonian();
    ~Hamiltonian();
    static double LocalEnergy(const arma::mat &, const int &, const int &, const double &);
    static double NumericalLocalEnergy(const arma::mat &, const int &, const int &, const double &, const double &, const double &);
    static double LocalEnergyInteraction(const arma::mat &, const int &, const int &, const double &, const double &, const double &);
    static double ParticleDistance(const arma::rowvec &, const arma::rowvec &);
    static arma::rowvec VectorSum(const arma::mat &, const int &, const int &, const double, const int &);
    static double DerivativeSum(const arma::mat &, const int &, const double, const int &);
    static double RepulsivePotential(const arma::mat &, const int &, const double &, const int &);
};

#endif // HAMILTONIAN_H

