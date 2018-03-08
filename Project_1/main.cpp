#define CATCH_CONFIG_RUNNER // Configure Catch to use this main, and not its own.
#include "inc/variationalmontecarlo.h"
#include "inc/catch.hpp"
#include "time.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>  // high resolution timing: http://en.cppreference.com/w/cpp/chrono/c/clock

#define TEST false // Change to true when testing and to false when running the program.


int RunCatchTests()
{
    return Catch::Session().run();
}

int main(int numberOfArguments, char *arguments[])
{
    if (TEST)
    {
        std::cout << "Running tests" << std::endl;
        return RunCatchTests();
    } else
    {
        // Default if there is no command line arguments
        int nParticles = 1;
        int nDimensions = 1;
        int nCycles = 1e6;
        double alpha = 0.5;
        double stepLength = 0.1;
        int cycleStepToFile = nCycles;

        int trials = 1;                                 // change to 10 when running timing

        // If command line arguments are defined
        if (numberOfArguments >= 2) { nParticles = std::atoi(arguments[1]); }
        if (numberOfArguments >= 3) { nDimensions = std::atoi(arguments[2]); }
        if (numberOfArguments >= 4) { nCycles = std::atoi(arguments[3]); }
        if (numberOfArguments >= 5) { alpha = std::atof(arguments[4]); }
        if (numberOfArguments >= 6) { stepLength = std::atof(arguments[5]); }
        if (numberOfArguments >= 7) { cycleStepToFile = std::atoi(arguments[6]); }
        if (numberOfArguments >= 8) { trials = std::atoi(arguments[7]); }

        std::cout << std::endl;
        std::cout << "Number of particles = " << nParticles << std::endl;
        std::cout << "Number of dimensions = " << nDimensions << std::endl;
        std::cout << "Number of cycles = " << nCycles << std::endl;
        std::cout << "Alpha = " << alpha << std::endl;
        std::cout << "Step length = " << stepLength << std::endl;
        std::cout << std::endl;

        VariationalMonteCarlo *VMC = new VariationalMonteCarlo();

        std::vector<double> timing = {};

        // TODO: Move timing to the MC cycles loop
        for (int i = 0; i < trials; i++)
        {
            // Start timing
            auto start_time = std::chrono::high_resolution_clock::now();

            VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, alpha, stepLength, cycleStepToFile);

            // Timing finished
            auto end_time = std::chrono::high_resolution_clock::now();

            timing.push_back(std::chrono::duration<double> (end_time - start_time).count());
        }
        double sum = 0;
        for (int i = 0; i < trials; i++)
        {
            sum += timing.at(i);
        }
        double averageTime = sum/timing.size();

        std::cout << std::endl;
        std::cout << "Average run time_chrono (Wall clock time): " << averageTime << " sec." << std::endl;
        std::cout << std::endl;

        return 0;
    }
}
