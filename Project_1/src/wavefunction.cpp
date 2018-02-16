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
        r[i] = sqrt(R[i][0]*R[i][0]+R[i][1]*R[i][1]+R[i][2]*R[i][2]);
        g[i] = exp(r[i]*r[i]);
        for (int j=0; j<nParticles; j++)
        {
            r[j] = sqrt(R[j][0]*R[j][0]+R[j][1]*R[j][1]+R[j][2]*R[j][2]);
            if (std::abs(r[i]-r[j])>a) {f[i][j] = 1-a/std::abs(r[i]-r[j]);}
            else {f[i][j] = 0;}
        }
    }
    //psi[i] = g*f;
}


Wavefunction::~Wavefunction()
{

}
