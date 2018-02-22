#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "armadillo"


class Wavefunction
{
public:
    Wavefunction(int nParticles, const arma::mat &);
    ~Wavefunction();
    //double Wavefunction(int nParticles, double** positionMatrix);
    double calculate_psi(int nParticles, const arma::mat &);
    double g(double position);
    double f(const arma::vec &position1, const arma::vec &position2, double a);
    double Hamiltonian(int nParticles, const arma::mat &);
    double localEnergy(int nParticles, const arma::mat &positionMatrix);



};

#endif // WAVEFUNCTION_H
