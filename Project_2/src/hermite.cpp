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
