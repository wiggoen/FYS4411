#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H


class Wavefunction
{
public:
    Wavefunction(int nParticles, double** positionMatrix);
    ~Wavefunction();
    //double Wavefunction(int nParticles, double** positionMatrix);
    double calculate_psi(int nParticles, double** positionMatrix);


};

#endif // WAVEFUNCTION_H
