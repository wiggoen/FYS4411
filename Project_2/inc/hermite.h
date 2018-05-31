#pragma once
#ifndef HERMITE_H
#define HERMITE_H


class Hermite
{
public:
    Hermite( void );
    ~Hermite( void );
    static double H(const int &n, const double &x);
    static double DerivativeHermite(const int &n, const double &x);
    static double DoubleDerivativeHermite(const int &n, const double &x);
};

#endif // HERMITE_H
