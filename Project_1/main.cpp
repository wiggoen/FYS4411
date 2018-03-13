#define CATCH_CONFIG_RUNNER // Configure Catch to use this main, and not its own.
#include "inc/variationalmontecarlo.h"
#include "inc/catch.hpp"


#define TEST false // Change to true when testing and to false when running the program.


int RunCatchTests()
{
    return Catch::Session().run();
}


int main(int numberOfArguments, char *arguments[])
{
    if (TEST)
    {
        std::cout << "Running tests..." << std::endl;
        return RunCatchTests();
    } else
    {
        // Default if there is no command line arguments
        int nParticles = 1;
        int nDimensions = 1;
        int nCycles = 1e6;
        double alpha = 0.5;
        double stepLength = 0.1;
        double dt = 0.01;               // Time step interval [0.001,0.01]
        int cycleStepToFile = nCycles;
        int trials = 1;                 // change to 10 when running timing

        // If command line arguments are defined
        if (numberOfArguments >= 2) { nParticles = std::atoi(arguments[1]); }
        if (numberOfArguments >= 3) { nDimensions = std::atoi(arguments[2]); }
        if (numberOfArguments >= 4) { nCycles = std::atoi(arguments[3]); }
        if (numberOfArguments >= 5) { alpha = std::atof(arguments[4]); }
        if (numberOfArguments >= 6) { stepLength = std::atof(arguments[5]); }
        if (numberOfArguments >= 7) { dt = std::atof(arguments[6]); }
        if (numberOfArguments >= 8) { cycleStepToFile = std::atoi(arguments[7]); }
        if (numberOfArguments >= 9) { trials = std::atoi(arguments[8]); }

        // Initialize VMC
        VariationalMonteCarlo *VMC = new VariationalMonteCarlo();


        // Setup for writing to file
        std::cout << std::endl;
        std::cout << "Particles " << " Dimensions " << "    Cycles " << " Alpha " << " Step_length "
                  << " Time_step " << " Time_[sec] " << " Energy " << " Energy_squared "
                  << " Variance " << " Acceptance_ratio " << std::endl;

        // Allocation
        arma::rowvec runVector;
        arma::mat runMatrix;

        // Run VMC
        for (int i = 0; i < trials; i++)
        {
            runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, alpha,
                                                      stepLength, dt, cycleStepToFile);
            runMatrix.insert_rows(i, runVector);
        }
        arma::rowvec columnSum;
        columnSum = arma::sum(runMatrix, 0);

        // Finding averages of trials
        double time = columnSum(0)/trials;
        double energy = columnSum(1)/trials;
        double energySquared = columnSum(2)/trials;
        double variance = columnSum(3)/trials;
        double acceptanceRatio = columnSum(4)/trials;

        std::cout << std::setw(9) << std::setprecision(3) << nParticles
                  << std::setw(12) << std::setprecision(3) << nDimensions
                  << std::setw(11) << std::setprecision(8) << nCycles
                  << std::setw(7) << std::setprecision(3) << alpha
                  << std::setw(13) << std::setprecision(3) << stepLength
                  << std::setw(11) << std::setprecision(6) << dt
                  << std::setw(12) << std::setprecision(6) << time
                  << std::setw(8) << std::setprecision(3) << energy
                  << std::setw(16) << std::setprecision(3) << energySquared
                  << std::setw(10) << std::setprecision(3) << variance
                  << std::setw(18) << std::setprecision(6) << acceptanceRatio << std::endl;

        return 0;
    }
}
