#include <iostream>
#include <fstream>
#include <string>
#include "inc/variationalmontecarlo.h"
#include "inc/json.hpp"

using json = nlohmann::json; // for convenience


int main()//int numberOfArguments, char *arguments[])
{
    /* The whole project path is needed to read json on mac */

    /* Who am I? */
    std::string Iam = "Trond";
    std::string projectFolder;
    if (Iam == "Trond") {
        std::cout << "Hi, Trond!" << std::endl;
        projectFolder = "/Users/trondwj/GitHub/FYS4411/Project_2/";
    } else if (Iam == "Line") {
        std::cout << "Hi, Line!" << std::endl;
        projectFolder = "/home/line/github/FYS4411/Project_2/";
    } else {
        std::cerr << "I don't know who you are.." << std::endl;
        exit(1);
    }


    /* READ PARAMETERS */
    std::ifstream inputFile(projectFolder+"parameters.json", std::ios::in);
    json parameter;
    inputFile >> parameter;
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
    std::string derivationType = parameter["Derivation type"];
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
                                UseJastrowFactor, UseImportanceSampling, UseFermionInteraction, derivationType, cycleType);
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


    /* Write to terminal */
    std::cout << std::endl;
    std::cout << "Jastrow_factor "  << " Importance_sampling " << " Fermion_interaction " << " MPI " << " Average_timing "
              << " Derivation_type " << "       Cycle_type " << std::endl;

    std::cout << std::setw(14) << Jastrow
              << std::setw(21) << Importance
              << std::setw(21) << Interaction
              << std::setw(5)  << Parallel
              << std::setw(16) << Timing
              << std::setw(17) << derivationType
              << std::setw(18) << cycleType << std::endl;


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
