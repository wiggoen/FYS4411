#include "inc/variationalmontecarlo.h"
#include "inc/wavefunction.h"
#include "inc/matrix.h"
#include <random>
#include <iomanip>
#include <iostream>

VariationalMonteCarlo::VariationalMonteCarlo(int rows, int columns)
{
    N = 1;
    Matrix matrix;

    Wavefunction wavefunc(N,matrix);

    int row = rows;     //number of particles
    int col = columns;  //dimensions
    int **matrix = makeMatrix(row, col);
    printMatrix(matrix, row, col);

    // Propose a new position R by moving one boson at the time
    // * Pick random boson
    int random_boson = rand() % row;
    // * Move boson
    for (int i = 0; i<col; i++)
    {
        matrix[random_boson][i] = random_double(1);
    }

    // Calculate new psi
    double psi = wavefunc(row,matrix);
    // Pick random number r in [0,1]
    double r = random_double(1);
    // Test if r is smaller or equal to |psi_T(R')|^2/|psi_T(R')|^2 ??
    double test_func = 1; //change this to be actual function above
    // If yes: accept new position
    if (r<test_func)
    {
        // Calculate delta E_L(R)
        // Update E = E + delta E_L
        // E^2 = E^2 + (delta E_L)^2
    }else
    {
        new_position = old_position;
    }
    return new_postition;

}

VariationalMonteCarlo::~VariationalMonteCarlo(int rows, int columns)
{

}

double random_double(fMax)
{
    double f = (double)rand() / RAND_MAX;
    return f*fMax;
}

int **VariationalMonteCarlo::makeMatrix(int rows, int columns)
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
void VariationalMonteCarlo::printMatrix(int *Matrix[], int rows, int columns)
{
    std::cout << std::endl;
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            std::cout << std::setw(10) << std::setprecision(3) << Matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
