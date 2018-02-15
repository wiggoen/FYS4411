#ifndef VARIATIONALMONTECARLO_H
#define VARIATIONALMONTECARLO_H


class variationalMonteCarlo
{
public:
    variationalMonteCarlo(int rows, int columns);
    ~variationalMonteCarlo();
    int **makeMatrix(int rows, int columns);
    void printMatrix(int *Matrix[], int rows, int columns);
private:
    int** Matrix;
    int rows;
    int columns;
};

#endif // VARIATIONALMONTECARLO_H
