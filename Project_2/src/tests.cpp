#ifndef CATCH_HPP
#include "inc/catch.hpp"
#include "inc/variationalmontecarlo.h"


VariationalMonteCarlo *VMC = new VariationalMonteCarlo();


TEST_CASE("Local energy two electrons", "[Hamiltonian, BF, IS]")
{
    int terminalizationFactor = 10000;
    int nParticles = 2;
    int nCycles = pow(2, 20) + terminalizationFactor;
    double alpha = 1.0;
    double beta = 0.0;
    double omega;
    double stepLength = 0.5;
    bool UseImportanceSampling;
    double timeStep = 0.001;
    bool UseJastrowFactor = false;
    bool UseFermionInteraction = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;

    arma::rowvec runVector;
    double energy;
    double variance;

    // Brute force
    UseImportanceSampling = false;
    omega = 1.0;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 2.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ) );


    omega = 0.5;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 1.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ) );


    omega = 0.1;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 0.2 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );


    omega = 0.05;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 0.1 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );

    omega = 0.01;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 0.02 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );


    // Importance sampling
    UseImportanceSampling = true;
    omega = 1.0;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 2.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ) );


    omega = 0.5;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 1.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ) );


    omega = 0.1;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 0.2 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );


    omega = 0.05;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 0.1 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );

    omega = 0.01;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 0.02 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );
}


TEST_CASE("Local energy six electrons", "[Hamiltonian, BF, IS]")
{
    int terminalizationFactor = 10000;
    int nParticles = 6;
    int nCycles = pow(2, 20) + terminalizationFactor;
    double alpha = 1.0;
    double beta = 0.0;
    double omega;
    double stepLength = 0.5;
    bool UseImportanceSampling;
    double timeStep = 0.001;
    bool UseJastrowFactor = false;
    bool UseFermionInteraction = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;


    arma::rowvec runVector;
    double energy;
    double variance;

    // Brute force
    UseImportanceSampling = false;
    omega = 1.0;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 10.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ) );


    // Importance sampling
    UseImportanceSampling = true;
    omega = 1.0;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 10.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ) );
}


TEST_CASE("Local energy twelve electrons", "[Hamiltonian, BF, IS and omega 1]")
{
    int terminalizationFactor = 10000;
    int nParticles = 12;
    int nCycles = pow(2, 20) + terminalizationFactor;
    double alpha = 1.0;
    double beta = 0.0;
    double omega;
    double stepLength = 0.5;
    bool UseImportanceSampling;
    double timeStep = 0.001;
    bool UseJastrowFactor = false;
    bool UseFermionInteraction = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;


    arma::rowvec runVector;
    double energy;
    double variance;

    // Brute force
    UseImportanceSampling = false;
    omega = 1.0;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 28.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );


    // Importance sampling
    UseImportanceSampling = true;
    omega = 1.0;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 28.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );
}


TEST_CASE("Local energy twenty electrons", "[Hamiltonian, BF, IS]")
{
    int terminalizationFactor = 10000;
    int nParticles = 20;
    int nCycles = pow(2, 20) + terminalizationFactor;
    double alpha = 1.0;
    double beta = 0.0;
    double omega;
    double stepLength = 0.5;
    bool UseImportanceSampling;
    double timeStep = 0.001;
    bool UseJastrowFactor = false;
    bool UseFermionInteraction = false;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;


    arma::rowvec runVector;
    double energy;
    double variance;

    // Brute force
    UseImportanceSampling = false;
    omega = 1.0;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 60.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-16) );


    // Importance sampling
    UseImportanceSampling = true;
    omega = 1.0;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy == Approx( 60.0 ) );
    variance = runVector(3);
    REQUIRE( variance == Approx( 0.0 ).margin(1e-13) );
}


TEST_CASE("Local energy two electrons, different alpha and beta", "[Hamiltonian, Jastrow, Interaction]")
{
    int terminalizationFactor = 10000;
    int nParticles = 2;
    int nCycles = pow(2, 20) + terminalizationFactor;
    double alpha;
    double beta;
    double omega = 1.0;
    double stepLength = 0.5;
    bool UseImportanceSampling;
    double timeStep = 0.001;
    bool UseJastrowFactor;
    bool UseFermionInteraction;
    bool UseAverageTiming = false;
    bool UseAnalyticalExpressions = true;
    bool UseNumericalPotentialEnergy = false;
    std::string cycleType = "MonteCarlo";
    int cycleStepToFile = 0;

    arma::rowvec runVector;
    double energy;
    double variance;

    alpha = 0.980576;
    beta = 0.434712;
    UseJastrowFactor = true;
    UseFermionInteraction = true;
    UseImportanceSampling = false;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy >= 3.0 );

    alpha = 0.774883;
    UseJastrowFactor = false;
    UseFermionInteraction = true;
    UseImportanceSampling = true;
    runVector = VMC->RunVMC(nParticles, nCycles, alpha, beta, omega, stepLength, timeStep, UseJastrowFactor,
                            UseImportanceSampling, UseFermionInteraction, UseAnalyticalExpressions,
                            UseNumericalPotentialEnergy, cycleType, cycleStepToFile, terminalizationFactor);
    energy = runVector(1);
    REQUIRE( energy > 3.0 );
}

#endif /* CATCH_HPP */
