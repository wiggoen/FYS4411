#ifndef CATCH_HPP
#include "inc/catch.hpp"
#include "inc/variationalmontecarlo.h"


VariationalMonteCarlo *VMC = new VariationalMonteCarlo();


TEST_CASE("Local energy two electrons", "[Hamiltonian]")
{
    int nParticles = 2;
    int nCycles = pow(2, 20);
    double alpha = 1;
    double beta = 0;
    double omega = 1;
    double spinParameter = 1;
    double stepLength = 1.5;
    bool UseImportanceSampling;
    double timeStep = 0;
    bool UseJastrowFactor = false;
    bool UseFermionInteraction = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;
    bool UseMPI = false;


    arma::rowvec runVector;
    double energy;
    double variance;

    UseImportanceSampling = false;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, spinParameter, stepLength, timeStep,
                            UseJastrowFactor, UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, UseMPI);
    energy = runVector(1);
    REQUIRE(energy == 2);
    variance = runVector(3);
    REQUIRE(variance == 0);

    UseImportanceSampling = true;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, spinParameter, stepLength, timeStep,
                            UseJastrowFactor, UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, UseMPI);
    energy = runVector(1);
    REQUIRE(energy == 2);
    variance = runVector(3);
    REQUIRE(variance == 0);
}

/*
TEST_CASE("Local energy six electrons", "[Hamiltonian]")
{
    int nParticles = 6;
    int nCycles = pow(2, 20);
    double alpha = 1;
    double beta = 0;
    double omega = 1;
    double spinParameter = 0;
    double stepLength = 0.5;
    bool UseImportanceSampling = false;
    double timeStep = 0;
    bool UseJastrowFactor = false;
    bool UseFermionInteraction = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;
    bool UseMPI = false;


    arma::rowvec runVector;
    double energy;
    double variance;

    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, spinParameter, stepLength, timeStep,
                            UseJastrowFactor, UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, UseMPI);
    energy = runVector(1);
    REQUIRE(energy == 10);
    variance = runVector(3);
    REQUIRE(variance == 0);
}


TEST_CASE("Local energy twelve electrons", "[Hamiltonian]")
{
    int nParticles = 12;
    int nCycles = pow(2, 20);
    double alpha = 1;
    double beta = 0;
    double omega = 1;
    double spinParameter = 0;
    double stepLength = 0.5;
    bool UseImportanceSampling = false;
    double timeStep = 0;
    bool UseJastrowFactor = false;
    bool UseFermionInteraction = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;
    bool UseMPI = false;


    arma::rowvec runVector;
    double energy;
    double variance;

    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, spinParameter, stepLength, timeStep,
                            UseJastrowFactor, UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, UseMPI);
    energy = runVector(1);
    REQUIRE(energy == 28);
    variance = runVector(3);
    REQUIRE(variance == 0);
}


TEST_CASE("Local energy twenty electrons", "[Hamiltonian]")
{
    int nParticles = 20;
    int nCycles = pow(2, 20);
    double alpha = 1;
    double beta = 0;
    double omega = 1;
    double spinParameter = 0;
    double stepLength = 0.5;
    bool UseImportanceSampling = false;
    double timeStep = 0;
    bool UseJastrowFactor = false;
    bool UseFermionInteraction = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;
    bool UseMPI = false;


    arma::rowvec runVector;
    double energy;
    double variance;

    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, spinParameter, stepLength, timeStep,
                            UseJastrowFactor, UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, UseMPI);
    energy = runVector(1);
    REQUIRE(energy == 60);
    variance = runVector(3);
    REQUIRE(variance == 0);
}
*/
#endif /* CATCH_HPP */
