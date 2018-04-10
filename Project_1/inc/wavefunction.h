#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>


class Wavefunction
{
public:
    Wavefunction();
    ~Wavefunction();
    static double TrialWaveFunction(const arma::mat &, const int &, const int &, const double &, const double &);
    static double TrialWaveFunctionInteraction(const arma::mat &, const int &, const int &, const double &, const double &, const double &);
    static void QuantumForce(const arma::mat &, arma::mat &, const double &);
    static void NumericalQuantumForce(const arma::mat &, arma::mat &, const int &, const int &, const double &, const double &, const double &);
    static void QuantumForceInteraction(const arma::mat &, arma::mat &, const double &, const int &, const int &, const double &, const int &);
    static double DerivativePsi(const arma::mat &, const int &, const int &, const double &);
};

#endif // WAVEFUNCTION_H
