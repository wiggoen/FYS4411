#pragma once
#ifndef HERMITE_H
#define HERMITE_H
#include <armadillo>


class Hermite
{
public:
    Hermite( void );
    ~Hermite( void );

    static double H(const int &n, const double &x);

    static double DerivativeHermite(const int &n, const double &x);

    static double DoubleDerivativeHermite(const int &n, const double &x);

    static double AlphaDerivativeHermite(const int &n, const double &x, const double &alpha, const double &omega);

    static arma::mat QuantumNumbers( void );

private:
};

#endif // HERMITE_H
