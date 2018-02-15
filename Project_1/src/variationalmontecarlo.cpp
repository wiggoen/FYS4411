#include "inc/variationalmontecarlo.h"
#include <random>
#include <iomanip>
#include <iostream>

variationalMonteCarlo::variationalMonteCarlo(int rows, int columns)
{

    int row = rows;
    int col = columns;
    int **matrix = makeMatrix(row, col);
    printMatrix(matrix, row, col);

    // Propose a new position R by moving one boson at the time
    // Calculate new psi
    // Pick random number r in [0,1]
    // Test if r is smaller or equal to |psi_T(R')|^2/|psi_T(R')|^2 ??
    // If yes: accept new position
        // Calculate delta E_L(R)
        // Update E = E + delta E_L
        // E^2 = E^2 + (delta E_L)^2

}

variationalMonteCarlo::~variationalMonteCarlo()
{

}

int **variationalMonteCarlo::makeMatrix(int rows, int columns)
{
    // Initialize the seed and call the Mersienne algorithm
    std::random_device rd;
    std::mt19937_64 gen(rd());

    // Set up the uniform distribution for x in [0, 1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);

    // Initialize matrix by dynamic memory allocation
    int **Matrix = new int *[rows];
    for (int i = 0; i < rows; i++)
    {
        Matrix[i] = new int[columns];
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
                Matrix[i][j] = -1;
            } else
            {
                Matrix[i][j] = randomNumber;
            }
        }
    }
    return Matrix;
}

// It is not recommended to print large matrices
// TODO: Set max value of print dimension
void variationalMonteCarlo::printMatrix(int *Matrix[], int rows, int columns)
{
    std::cout << std::endl;
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            std::cout << std::setw(3) << Matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
