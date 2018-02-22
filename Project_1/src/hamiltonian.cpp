#include "inc/hamiltonian.h"
#include "armadillo"

Hamiltonian::Hamiltonian(int nParticles, const arma::mat &)
{
    // H = sum(i-N){-hbar^2/2m del_i^2 + V_ext(r_i)}
    // + sum(i<j-N){V_int(r_i,r_j)}

    // For now this only calculates the spherical potential
    double sum1, sum2;
    int N = nParticles;
    for (int i=0; i<N; i++)
    {
        //wtf, don't know how to implement this...
    }
}

Hamiltonian::~Hamiltonian()
{

}
