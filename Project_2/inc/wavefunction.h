#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include <armadillo>


class Wavefunction
{
public:
    Wavefunction( void );
    ~Wavefunction( void );
    static double TrialWaveFunction(const arma::mat &r, const double &alpha, const double &beta,
                                    const double &omega, const double &a, const double &constant);
    static void QuantumForce(const arma::mat &r, arma::mat &QForce, const double &alpha, const double beta, const double a, const double omega);
    /*
    void NumericalQuantumForce(const arma::mat &r, arma::mat &QForce, const int &nParticles,
                                    const int &nDimensions, const double &alpha, const double &stepLength,
                                    const double &beta, const double &omega, const double &a, const double &constant);*/
};


#endif /* WAVEFUNCTION_H */
