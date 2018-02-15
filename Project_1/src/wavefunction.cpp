#include "inc/wavefunction.h"
#include <math.h>

Wavefunction::Wavefunction(int nParticles, double** positionMatrix)
{
    //psi = prod(i){g(alpha,beta,r_i)}
    // * prod(i<j){f(a,|r_i-r_j|}

    double* psi = new double(100);



    double g; // = exp(x^2+y^2+z^2);
    double f;

    int i = 0;
    double* r = new double(nParticles);

    for (int i=0; i<nParticles; i++)
    {
        r[i] = sqrt(positionMatrix[i,0],positionMatrix[i,1],positionMatrix[i,3]);

        g = exp(r[i]*r[i]);
        for (int j=0; j<nParticles; j++)
        {
            r[i] = sqrt(positionMatrix[j,0],positionMatrix[j,1],positionMatrix[j,3]);
                        if (abs(r[i]-r[j])>a) {f = 1-a/abs(r[i]-r[j]);}
                        else {f = 0;}
            }
        }
    }


Wavefunction::~Wavefunction()
{

}
