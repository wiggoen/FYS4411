#include "inc/hermite.h"
#include <iostream>


Hermite::Hermite( void )
{

}


Hermite::~Hermite( void )
{

}


double Hermite::H(const int &n, const double &x)
{
    switch(n)
    {
    case 0: return 1;
    case 1: return 2*x;
    case 2: return 4*x*x - 2;
    case 3: return 8*x*x*x - 12*x;
    case 4: return 16*x*x*x*x - 48*x*x + 12;
    default: std::cerr << "Something went wrong in the Hermite function" << std::endl; exit(1);
    }
}


double Hermite::DerivativeHermite(const int &n, const double &x)
{
    switch(n)
    {
    case 0: return 0;
    case 1: return 2;
    case 2: return 8*x;
    case 3: return 24*x*x - 12;
    case 4: return 64*x*x*x - 96*x;
    default: std::cerr << "Something went wrong in the DerivateHermite function" << std::endl; exit(1);
    }
}


double Hermite::DoubleDerivativeHermite(const int &n, const double &x)
{
    switch(n)
    {
    case 0: return 0;
    case 1: return 0;
    case 2: return 8;
    case 3: return 48*x;
    case 4: return 192*x*x - 96;
    default: std::cerr << "Something went wrong in the DoubleDerivateHermite function" << std::endl; exit(1);
    }
}

double Hermite::AlphaDerivativeHermite(const int &n, const double &x, const double &alpha, const double &omega)
{
    double sqrtAlphaOmega     = sqrt(alpha*omega);
    switch(n)
    {
    case 0: return 0;
    case 1: return omega*x/sqrtAlphaOmega;
    case 2: return 4*omega*x*x;
    case 3: return 12*sqrtAlphaOmega*omega*x*x*x;
    case 4: return 32*alpha*omega*omega*x*x*x*x - 48*omega*x*x;
    default: std::cerr << "Something went wrong in the AlphaDerivativeHermite function" << std::endl; exit(1);
    }
}


arma::mat Hermite::QuantumNumbers( void ) {
    /* Quantum numbers for up to 20 electrons */
    const arma::mat QNumbers = { {0, 0},
                                 {1, 0},
                                 {0, 1},
                                 {2, 0},
                                 {1, 1},
                                 {0, 2},
                                 {3, 0},
                                 {2, 1},
                                 {1, 2},
                                 {0, 3} };
    return QNumbers;
}
