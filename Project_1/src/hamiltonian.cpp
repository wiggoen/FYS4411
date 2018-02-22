#include "inc/hamiltonian.h"
#include "armadillo"
#include "inc/wavefunction.h"




Hamiltonian::Hamiltonian(int nParticles, const arma::vec &x)
{
    // H = sum(i-N){-hbar^2/2m del_i^2 + V_ext(r_i)}
    // + sum(i<j-N){V_int(r_i,r_j)}


    wf = new Wavefunction(0,x);

    // For now this only calculates the spherical potential
    double sum1, sum2;
    double h = 0.0001;
    int N = nParticles;
    double alpha = 0.5;
    double V_ext = 0;
    for (int i=0; i<N; i++)
    {
        sum1 += - alpha * DoubleDerivative(N,h,x) + V_ext;
    }
}

double Hamiltonian::Derivative(int nParticles, double h, arma::vec x)
{
    double dude;
    for (int i=0; i<nParticles; i++)
    {
        dude = (wf->g(x(i)+h)-wf->g(x(i)-h))/(2*h);
    }
    return dude;
}

double Hamiltonian::DoubleDerivative(int nParticles, double h, const arma::vec x)
{
    double dude;
    for (int i=0; i<nParticles; i++)
    {
        dude = (wf->g(x(i)+h) + wf->g(x(i)-h) -2*wf->g(x(i)))/(h*h);
    }
    return dude;
}


Hamiltonian::~Hamiltonian()
{

}
