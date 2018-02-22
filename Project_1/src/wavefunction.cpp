#include "inc/wavefunction.h"
#include "inc/matrix.h"
#include <math.h>
#include <cmath>
#include "armadillo"

Wavefunction::Wavefunction(int nParticles, const arma::mat &)
{
    //Matrix matrix;
    double calculate_psi(int nParticles, double** positionMatrix);
    double g(double* position);
    double f(double* position1, double* position2, double a);
}



double Wavefunction::g(const arma::vec &position)
{
    arma::vec p = position;
    double r = sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2));
    return exp(r*r);
}

double Wavefunction::f(const arma::vec &position1, const arma::vec &position2, double a)
{
    arma::vec p1 = position1;
    arma::vec p2 = position2;
    double r1 = sqrt(p1(0)*p1(0)+p1(1)*p1(1)+p1(2)*p1(2));
    double r2 = sqrt(p2(0)*p2(0)+p2(1)*p2(1)+p2(2)*p2(2));
    if (std::abs(r1-r2) > a)
    {
        return 1-a/(std::abs(r1-r2));
    } else {return 0;}
}

double Wavefunction::calculate_psi(int nParticles, const arma::mat &positionMatrix)
{
    //Matrix matrix;
    arma::mat R = positionMatrix;
    double a = 1;
    double psi;

    for (int i=0; i<nParticles; i++)
    {
        double g_i = g(R.row(i));
        for (int j=0; j<nParticles; j++)
        {
            double f_ij = f(R.row(i),R.row(j),a);

            psi *= g_i*f_ij;
        }
    }
    return psi;
}



double Wavefunction::localEnergy(int nParticles, const arma::mat &positionMatrix)
{
    // (1/psi)*H*psi
    int N = nParticles;
    double psi = calculate_psi(N,positionMatrix);

}


Wavefunction::~Wavefunction()
{

}
