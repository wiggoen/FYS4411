#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H


class VariationalMonteCarlo
{

private:
    int** Matrix = new int*[10];

public:
    VariationalMonteCarlo(int rows, int columns, double** positionMatrix);
    ~VariationalMonteCarlo();
    static void printMatrix(int *Matrix[], int rows, int columns);
    int **makeMatrix(int rows, int columns);// {return Matrix;}

};

#endif // VARIATIONALMONTECARLO_H
