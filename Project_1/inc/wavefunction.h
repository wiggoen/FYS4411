#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H


class Wavefunction
{
public:
    Wavefunction(int nParticles, double** positionMatrix);
    ~Wavefunction();
    //double Wavefunction(int nParticles, double** positionMatrix);
    double calculate_psi(int nParticles, double** positionMatrix);
    double g(double* position);
    double f(double* position1, double* position2, double a);
    double Hamiltonian(int nParticles, double** positionMatrix);
    double localEnergy(int nParticles, double** positionMatrix);



};

#endif // WAVEFUNCTION_H
