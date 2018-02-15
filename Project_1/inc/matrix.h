#ifndef MATRIX_H
#define MATRIX_H


class Matrix
{
public:
    Matrix(int rows, int columns);
    ~Matrix();
    int **makeMatrix(int rows, int columns);
    void printMatrix(int *matrix[], int rows, int columns);
private:
    int** matrix;
    int rows;
    int columns;
};

#endif // MATRIX_H
