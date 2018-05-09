#include <iostream>
#include <fstream>
#include <string>
#include "inc/variationalmontecarlo.h"
#include "inc/json.hpp"

using json = nlohmann::json; // for convenience


int main(int numberOfArguments, char *arguments[])
{
    /* The whole project path is needed to read json on mac */

    /* Who am I? */
    std::string Iam = "Line";
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
    int nDimensions = parameter["nDimensions"];
    int nCycles = parameter["nCycles"];
    double alpha = parameter["alpha"];
    double beta = parameter["beta"];
    double omega = parameter["omega"];
    double a = parameter["a"];
    double stepLength = parameter["stepLength"];
    double constant = parameter["constant"];


    /* Initialize VMC */
    VariationalMonteCarlo *VMC = new VariationalMonteCarlo();

    /* Allocation of information to be printed */
    arma::rowvec runVector;
    arma::mat runMatrix;

    /* Run VMC */
    int trials = 1;                     // For timing purposes
    for (int i = 0; i < trials; i++)
    {
        runVector = VMC->RunMonteCarloIntegration(nParticles, nCycles, alpha, beta, omega, a, stepLength, constant);
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
    std::cout << std::endl;

    std::cout << "Particles "  << " Dimensions " << "    Cycles " << " Alpha " << " Step_length "
              << " Time_step " << " Time_[sec] " << "   Energy "  << " Energy_squared "
              << " Variance "  << " Acceptance_ratio " << std::endl;

    /* Write to terminal */
    std::cout << std::setw(9)  << std::setprecision(3) << nParticles
              << std::setw(12) << std::setprecision(3) << nDimensions
              << std::setw(11) << std::setprecision(8) << nCycles
              << std::setw(7)  << std::setprecision(3) << alpha
              << std::setw(13) << std::setprecision(3) << stepLength
              << std::setw(11) << std::setprecision(6) << "Ignore" //timeStep
              << std::setw(12) << std::setprecision(6) << runTime
              << std::setw(10) << std::setprecision(3) << energy
              << std::setw(16) << std::setprecision(3) << energySquared
              << std::setw(10) << std::setprecision(3) << variance
              << std::setw(18) << std::setprecision(6) << acceptanceRatio
              << std::endl;

    return 0;
}
