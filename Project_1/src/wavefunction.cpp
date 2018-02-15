#include "inc/wavefunction.h"
#include <math.h>
#include <cmath>

Wavefunction::Wavefunction(int nParticles, double** positionMatrix)
{
    //psi = prod(i){g(alpha,beta,r_i)}
    // * prod(i<j){f(a,|r_i-r_j|}

    double* psi = new double(100);
    double** R = positionMatrix;



    double g; // = exp(x^2+y^2+z^2);
    double f;

    double a = 1;
    int i = 0;
    double* r = new double(nParticles);

    for (int i=0; i<nParticles; i++)
    {
        r[i] = sqrt(R[i][0]*R[i][0]+R[i][1]*R[i][1]+R[i][2]*R[i][2]);

        g = exp(r[i]*r[i]);
        for (int j=0; j<nParticles; j++)
        {
            r[j] = sqrt(R[j][0]*R[j][0]+R[j][1]*R[j][1]+R[j][2]*R[j][2]);
                        if (std::abs(r[i]-r[j])>a) {f = 1-a/std::abs(r[i]-r[j]);}
                        else {f = 0;}
            }
        }
    }


Wavefunction::~Wavefunction()
{

}
