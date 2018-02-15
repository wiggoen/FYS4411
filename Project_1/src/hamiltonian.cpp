#include "inc/hamiltonian.h"

Hamiltonian::Hamiltonian()
{
    // H = sum(i-N){-hbar^2/2m del_i^2 + V_ext(r_i)}
    // + sum(i<j-N){V_int(r_i,r_j)}
}

Hamiltonian::~Hamiltonian()
{

}
