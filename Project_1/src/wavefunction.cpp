#include "inc/wavefunction.h"
#include "inc/matrix.h"
#include <math.h>
#include <cmath>

Wavefunction::Wavefunction(int nParticles, double** positionMatrix)
{
    //psi = prod(i){g(alpha,beta,r_i)}
    // * prod(i<j){f(a,|r_i-r_j|}

    Matrix matrix;

    int N = 100;
    double* psi = new double(N);

    double** R = positionMatrix;
    double* g = new double(N);
    double** f = matrix.makeMatrix(nParticles, 2);

    double a = 1;
    int i = 0;
    double* r = new double(nParticles);

    for (int i=0; i<nParticles; i++)
    {
        for (int j=0; j<nParticles; j++)
        {

        }
    }
    //psi[i] = g*f;
}

double g(double* position)
{
    double* p = position;
    double r = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    return exp(r*r);
}

double f(double* position1, double* position2, double a)
{
    double* p1 = position1;
    double* p2 = position2;
    double r1 = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
    double r2 = sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);
    if (std::abs(r1-r2) > a)
    {
        return 1-a/(std::abs(r1-r2));
    } else {return 0;}
}


Wavefunction::~Wavefunction()
{

}
