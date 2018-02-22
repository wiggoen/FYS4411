#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H
#include "armadillo"


class VariationalMonteCarlo
{

private:
    int** Matrix = new int*[10];

public:
    VariationalMonteCarlo(int rows, int columns, const arma::mat &);
    ~VariationalMonteCarlo();
    //static void printMatrix(int *Matrix[], int rows, int columns);
    int **makeMatrix(int rows, int columns);// {return Matrix;}
    static void vmc(int rows, int columns, const arma::mat &);


};

#endif // VARIATIONALMONTECARLO_H
