#define CATCH_CONFIG_RUNNER // Configure Catch to use this main, and not its own.
#include "inc/variationalmontecarlo.h"
#include "inc/catch.hpp"
#include "time.h"
#include <iostream>
#include <vector>
#include <chrono>  // high resolution timing: http://en.cppreference.com/w/cpp/chrono/c/clock

#define TEST false // Change to true when testing and to false when running the program.


int runCatchTests()
{
    return Catch::Session().run();
}

int main()
{
    if (TEST)
    {
        std::cout << "Running tests" << std::endl;
        return runCatchTests();
    } else
    {
        int nDimensions = 1;
        int nParticles = 1;
        int nCycles = 1e6;
        double alpha = 0.5;
        double stepLength = 0.1;

        VariationalMonteCarlo *solver = new VariationalMonteCarlo();

        std::vector<double> timing = {};
        std::vector<double> timing_chrono = {};
        int trials = 1;

        // Timing the algorithm
        clock_t start, finish;

        for (int i = 0; i < trials; i++)
        {
            // Start timing
            start = clock();
            auto start_time = std::chrono::high_resolution_clock::now();

            solver->runMonteCarloIntegration(nParticles, nDimensions, nCycles, alpha, stepLength);

            // Timing finished
            finish = clock();
            auto end_time = std::chrono::high_resolution_clock::now();

            timing.push_back(double (finish - start)/CLOCKS_PER_SEC);
            timing_chrono.push_back(std::chrono::duration<double> (end_time - start_time).count());
        }
        double sum = 0;
        double sum_chrono = 0;
        for (int i = 0; i < trials; i++)
        {
            sum += timing.at(i);
            sum_chrono += timing_chrono.at(i);
        }
        double averageTime = sum/timing.size();
        double averageTime_chrono = sum_chrono/timing_chrono.size();

        std::cout << std::endl << "Average run time (CPU time): " << averageTime << " sec." << std::endl;
        std::cout << "Average run time_chrono (Wall clock time): " << averageTime_chrono << " sec." << std::endl;

        return 0;
    }
}
