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
        int nParticles      = 6;
        int nDimensions     = 3;
        int nCycles         = 1e6;//nParticles;  // To make a total of 1 million samples
        double alpha        = 0.5;
        double stepLength   = 0.1;             //                         << Use 0.1 or 0.2
        double timeStep     = 0.01;            // Interval [0.001,0.01]
        int cycleStepToFile = nCycles;
        int trials          = 1;               // change to 10 when running timing

        // CHOOSE SAMPLING METHOD                       <<< --- CHOOSE ONLY ONE
        std::string samplingType = "BruteForce";
        //std::string samplingType = "Importance";

        // CHOOSE INTEGRATION METHOD                    <<< --- CHOOSE ONLY ONE
        //std::string derivationType = "Analytical";
        //std::string derivationType = "Numerical";
        std::string derivationType = "Interaction";

        // CHOOSE CYCLE TYPE                            <<< --- CHOOSE ONLY ONE
        std::string cycleType = "MonteCarlo";
        //std::string cycleType = "SteepestDescent";

        // If command line arguments are defined
        if (numberOfArguments >= 2)  { nParticles      = std::atoi(arguments[1]); }
        if (numberOfArguments >= 3)  { nDimensions     = std::atoi(arguments[2]); }
        if (numberOfArguments >= 4)  { nCycles         = std::atoi(arguments[3]); }
        if (numberOfArguments >= 5)  { alpha           = std::atof(arguments[4]); }
        if (numberOfArguments >= 6)  { stepLength      = std::atof(arguments[5]); }
        if (numberOfArguments >= 7)  { timeStep        = std::atof(arguments[6]); }
        if (numberOfArguments >= 8)  { cycleStepToFile = std::atoi(arguments[7]); }
        if (numberOfArguments >= 9)  { trials          = std::atoi(arguments[8]); }
        if (numberOfArguments >= 10) {
            samplingType = std::string(arguments[9]);
            // Renaming sampling type to run with the program implementation
            if      (samplingType == "bf")  { samplingType = "BruteForce"; }
            else if (samplingType == "im")  { samplingType = "Importance"; }
            else {
                std::cerr << "Error: You have to specify sampling type. Sampling type '" << samplingType << "' is not valid." << "\n"
                          << "Options are 'bf' (brute force) or 'im' (importance sampling)." << std::endl;
                exit(1);
            }
        }
        if (numberOfArguments >= 11) {
            derivationType = std::string(arguments[10]);
            // Renaming integration type to run with the program implementation
            if      (derivationType == "ana")  { derivationType = "Analytical"; }
            else if (derivationType == "num")  { derivationType = "Numerical"; }
            else if (derivationType == "int")  { derivationType = "Interaction"; }
            else {
                std::cerr << "Error: You have to specify integration type. Integration type '" << derivationType << "' is not valid." << "\n"
                          << "Options are 'ana' (analytic) or 'num' (numerical) or 'int' (interaction)." << std::endl;
                exit(1);
            }
        }

        // Initialize VMC
        VariationalMonteCarlo *VMC = new VariationalMonteCarlo();

        // Allocation
        arma::rowvec runVector;
        arma::mat runMatrix;

        // Run VMC
        for (int i = 0; i < trials; i++)
        {
            runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, alpha,
                                                      stepLength, timeStep, cycleStepToFile,
                                                      samplingType, derivationType, cycleType);
            runMatrix.insert_rows(i, runVector);
        }

        // TODO: IS THIS IN USE? OR CAN IT BE REMOVED?
        // Run steepest descend to find best alpha
        //double bestAlpha = VMC->SteepestDescent(nParticles, nDimensions);
        //runMatrix.insert_rows(i, runVector);
        //std::cout << "Steepest descend yields best alpha: alpha = " << bestAlpha << std::endl;

        arma::rowvec columnSum;
        columnSum = arma::sum(runMatrix, 0);

        // Finding averages of trials
        double runTime = columnSum(0)/trials;
        double energy = columnSum(1)/trials;
        double energySquared = columnSum(2)/trials;
        double variance = columnSum(3)/trials;
        double acceptanceRatio = columnSum(4)/trials;

        // Setup for writing to file
        std::cout << std::endl;
        std::cout << "Particles " << " Dimensions " << "    Cycles " << " Alpha " << " Step_length "
                  << " Time_step " << " Time_[sec] " << "   Energy " << " Energy_squared "
                  << " Variance " << " Acceptance_ratio " << std::endl;

        std::cout << std::setw(9)  << std::setprecision(3) << nParticles
                  << std::setw(12) << std::setprecision(3) << nDimensions
                  << std::setw(11) << std::setprecision(8) << nCycles
                  << std::setw(7)  << std::setprecision(3) << alpha
                  << std::setw(13) << std::setprecision(3) << stepLength
                  << std::setw(11) << std::setprecision(6) << timeStep
                  << std::setw(12) << std::setprecision(6) << runTime
                  << std::setw(10) << std::setprecision(3) << energy
                  << std::setw(16) << std::setprecision(3) << energySquared
                  << std::setw(10) << std::setprecision(3) << variance
                  << std::setw(18) << std::setprecision(6) << acceptanceRatio << std::endl;

        return 0;
    }
}
