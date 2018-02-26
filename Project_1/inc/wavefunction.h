#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H


class Wavefunction
{
public:
    Wavefunction(const arma::mat &r);
    ~Wavefunction();
    double waveFunction(const arma::mat &r);


};

#endif // WAVEFUNCTION_H
