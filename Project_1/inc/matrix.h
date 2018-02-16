#ifndef MATRIX_H
#define MATRIX_H


class Matrix
{
public:
    Matrix();
    ~Matrix();
    double **makeMatrix(int rows, int columns);
    void printMatrix(double *matrix[]);
    double **matrix;
private:
    int rows;
    int columns;
};

#endif // MATRIX_H
