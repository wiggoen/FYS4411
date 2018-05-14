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
    double spinParameter = 0;
    double stepLength = 1.5;
    double timeStep = 0;
    bool UseJastrowFactor = false;
    bool UseImportanceSampling = false;
    bool UseFermionInteraction = false;
    bool UseMPI = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";


    arma::rowvec runVector;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, spinParameter, stepLength, timeStep,
                            UseJastrowFactor, UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType);
    double energy = runVector(1);
    REQUIRE(energy == 2);
    double variance = runVector(3);
    REQUIRE(variance == 0);
}

#endif /* CATCH_HPP */
