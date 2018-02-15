#ifndef MATRIX_H
#define MATRIX_H


class Matrix
{
public:
    Matrix(int rows, int columns);
    ~Matrix();
    double **makeMatrix(int rows, int columns);
    void printMatrix(double *matrix[], int rows, int columns);
    double **matrix;
private:
    int rows;
    int columns;
};

#endif // MATRIX_H
