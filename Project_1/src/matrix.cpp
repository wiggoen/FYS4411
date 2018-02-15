#include "inc/matrix.h"
#include <random>
#include <iomanip>
#include <iostream>

Matrix::Matrix(int rows, int columns)
{
    int row = rows;
    int col = columns;
    int **matrix = makeMatrix(row, col);
    printMatrix(matrix, row, col);
}

Matrix::~Matrix()
{

}

int **Matrix::makeMatrix(int rows, int columns)
{
    // Initialize the seed and call the Mersienne algorithm
    std::random_device rd;
    std::mt19937_64 gen(rd());

    // Set up the uniform distribution for x in [0, 1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);

    // Initialize matrix by dynamic memory allocation
    int **matrix = new int *[rows];
    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new int[columns];
    }
    /* Deallocation of matrix
     * for (int i = 0; i < rows; i++) {
     *     delete [] Matrix[i];
     * delete [] Matrix;
     * }
    */

    // Random numbers or ground state in matrix
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < columns; j++) {
            int randomNumber = (int) (RandomNumberGenerator(gen)*2);
            std::cout << "randomNumber = " << randomNumber << std::endl;
            if (randomNumber == 0)
            {
                matrix[i][j] = -1;
            } else
            {
                matrix[i][j] = randomNumber;
            }
        }
    }
    return matrix;
}

// It is not recommended to print large matrices
// TODO: Set max value of print dimension
void Matrix::printMatrix(int *matrix[], int rows, int columns)
{
    std::cout << std::endl;
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            std::cout << std::setw(3) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
