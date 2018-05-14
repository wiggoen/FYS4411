#define CATCH_CONFIG_RUNNER /* Configure Catch to use this main, and not its own */
#include "inc/variationalmontecarlo.h"
#include "inc/catch.hpp"
#include "inc/json.hpp"


using json = nlohmann::json; // for convenience

int RunCatchTests()
{
    return Catch::Session().run();
}


int main(int argumentCount, char *argumentVector[])
{
    /* Allocating input file for parameters */
    std::ifstream inputFile;

    /* Check the number of command line arguments */
    if (argumentCount < 2) {
        /* Who am I? */
        std::string Iam = "Trond";

        /* The whole project path is needed to read json file */
        std::string projectFolder;
        if (Iam == "Trond") {
            std::cout << "Hi, Trond!" << std::endl;
            projectFolder = "/Users/trondwj/GitHub/FYS4411/Project_2/";
            inputFile.open(projectFolder+"parameters.json", std::ios::in);
        } else if (Iam == "Line") {
            std::cout << "Hi, Line!" << std::endl;
            projectFolder = "/home/line/github/FYS4411/Project_2/";
            inputFile.open(projectFolder+"parameters.json", std::ios::in);
        } else {
            std::cerr << "You have to give the parameter file as command line argument." << std::endl;
            std::cerr << "Compile: g++ -std=c++11 -O3 -L '/armadillo/lib/' -I '/armadillo/include/' -I 'project folder' -o main.out main.cpp src/*.cpp -larmadillo" << std::endl;
            std::cerr << "Remember to include the paths of armadillo library and headers as stated above, and in addition the project folder." << std::endl;
            std::cerr << "Usage: ./main.out parameters.json" << std::endl;
            exit(1);
        }
    } else
    {
        /* Input file from command line arguments */
        inputFile.open(argumentVector[1], std::ios::in);
    }

    /* READ PARAMETERS */
    json parameter;
    inputFile >> parameter;
    bool TEST = parameter["Run tests"];
    if (TEST) {
        std::cout << "Running tests..." << std::endl;
        return RunCatchTests();
    } else
    {
        int nParticles = parameter["nParticles"];
        int nCycles = parameter["nCycles"];
        double alpha = parameter["alpha"];
        double beta = parameter["beta"];
        double omega = parameter["omega"];
        double spinParameter = parameter["spin parameter"];
        double stepLength = parameter["stepLength"];
        double timeStep = parameter["timeStep"];
        bool UseJastrowFactor = parameter["Use Jastrow factor"];
        bool UseImportanceSampling = parameter["Use Importance sampling"];
        bool UseFermionInteraction = parameter["Use Fermion interaction"];
        bool UseMPI = parameter["Use MPI"];
        bool UseAverageTiming = parameter["Use Average timing"];
        bool UseAnalyticalExpressions = parameter["Use analytical expressions"];
        bool UseNumericalPotentialEnergy = parameter["Use numerical potential energy"];
        std::string cycleType = parameter["Cycle type"];


        /* Initialize VMC */
        VariationalMonteCarlo *VMC = new VariationalMonteCarlo();

        /* Allocation of information to be printed */
        arma::rowvec runVector;
        arma::mat runMatrix;

        /* For timing purposes */
        int trials;
        if (!UseAverageTiming)  /* No average timing */
        {
            trials = 1;
        } else
        {
            trials = 10;
        }

        /* Run VMC */
        for (int i = 0; i < trials; i++)
        {
            runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, spinParameter, stepLength, timeStep,
                                    UseJastrowFactor, UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                                    UseNumericalPotentialEnergy, cycleType);
            runMatrix.insert_rows(i, runVector);
        }

        /* Used for running several trials */
        arma::rowvec columnSum;
        columnSum = arma::sum(runMatrix, 0);

        /* Finding averages of trials */
        double runTime = columnSum(0)/trials;
        double energy = columnSum(1)/trials;
        double energySquared = columnSum(2)/trials;
        double variance = columnSum(3)/trials;
        double acceptanceRatio = columnSum(4)/trials;

        /* Setup for printing to terminal */
        int nDimensions = 2;

        std::string Jastrow;
        UseJastrowFactor ? Jastrow = "On" : Jastrow = "Off";  /* (condition) ? (if_true) : (if_false) */

        std::string Importance;
        UseImportanceSampling ? Importance = "On" : Importance = "Off";

        std::string Interaction;
        UseFermionInteraction ? Interaction = "On" : Interaction = "Off";

        std::string Parallel;
        UseMPI ? Parallel = "On" : Parallel = "Off";

        std::string Timing;
        UseAverageTiming ? Timing = "On" : Timing = "Off";

        std::string AnalyticExpressions;
        UseAnalyticalExpressions ? AnalyticExpressions = "On" : AnalyticExpressions = "Off";

        std::string NumericalPotentialEnergy;
        UseNumericalPotentialEnergy ? NumericalPotentialEnergy = "On" : NumericalPotentialEnergy = "Off";

        /* Write to terminal */
        std::cout << std::endl;
        std::cout << "Jastrow_factor "  << " Importance_sampling " << " Fermion_interaction " << " MPI " << " Average_timing "
                  << " Analytical_expressions " << " Numerical_potential_energy " << "      Cycle_type " << std::endl;

        std::cout << std::setw(14) << Jastrow
                  << std::setw(21) << Importance
                  << std::setw(21) << Interaction
                  << std::setw(5)  << Parallel
                  << std::setw(16) << Timing
                  << std::setw(24) << AnalyticExpressions
                  << std::setw(28) << NumericalPotentialEnergy
                  << std::setw(17) << cycleType << std::endl;


        std::cout << "Particles "  << " Dimensions " << "    Cycles " << " Alpha " << " Beta " << " Omega "
                  << " Spin_parameter "  << " Step_length " << " Time_step " << " Time_[sec] " << "      Energy "
                  << " Energy_squared " << " Variance "  << " Acceptance_ratio " << std::endl;

        std::cout << std::setw(9)  << std::setprecision(3) << nParticles
                  << std::setw(12) << std::setprecision(3) << nDimensions
                  << std::setw(11) << std::setprecision(8) << nCycles
                  << std::setw(7)  << std::setprecision(3) << alpha
                  << std::setw(6)  << std::setprecision(3) << beta
                  << std::setw(7)  << std::setprecision(3) << omega
                  << std::setw(16) << std::setprecision(3) << spinParameter
                  << std::setw(13) << std::setprecision(3) << stepLength
                  << std::setw(11) << std::setprecision(6) << timeStep
                  << std::setw(12) << std::setprecision(6) << runTime
                  << std::setw(13) << std::setprecision(6) << energy
                  << std::setw(16) << std::setprecision(6) << energySquared
                  << std::setw(10) << std::setprecision(3) << variance
                  << std::setw(18) << std::setprecision(6) << acceptanceRatio << std::endl;

        return 0;
    }
}
