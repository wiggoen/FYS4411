#ifndef CATCH_HPP
#include "inc/catch.hpp"
#include "inc/variationalmontecarlo.h"


VariationalMonteCarlo *VMC = new VariationalMonteCarlo();

TEST_CASE("Local energy, alpha 0.5", "[Hamiltonian]")
{
    int nCycles = 1e6;
    double alpha = 0.5;
    double stepLength = 0.1;
    double dt = 0.01;
    int cycleStepToFile = 1e6;
    double energy = 0;
    double variance = 0;
    double zero = 0;

    arma::rowvec runVector;
    runVector = VMC->RunMonteCarloIntegration(1, 1, nCycles, alpha, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy == 0.5);
    variance = runVector(3);
    REQUIRE(variance == zero);

    runVector = VMC->RunMonteCarloIntegration(1, 2, nCycles, alpha, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy == 1.0);
    variance = runVector(3);
    REQUIRE(variance == zero);

    runVector = VMC->RunMonteCarloIntegration(1, 3, nCycles, alpha, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy == 1.5);
    variance = runVector(3);
    REQUIRE(variance == zero);

    runVector = VMC->RunMonteCarloIntegration(10, 1, nCycles, alpha, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy == 5);
    variance = runVector(3);
    REQUIRE(variance == zero);

    runVector = VMC->RunMonteCarloIntegration(10, 2, nCycles, alpha, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy == 10);
    variance = runVector(3);
    REQUIRE(variance == zero);

    runVector = VMC->RunMonteCarloIntegration(10, 3, nCycles, alpha, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy == 15);
    variance = runVector(3);
    REQUIRE(variance == zero);
}


TEST_CASE("Local energy, 1 particle, 1D, alpha != 0.5", "[Hamiltonian, alpha not equal to 0.5, 1 particle, 1 dimension]")
{
    int nParticles = 1;
    int nDimensions = 1;
    int nCycles = 1e6;
    double stepLength = 0.1;
    double dt = 0.01;
    int cycleStepToFile = 1e6;

    double energy = 0;
    double localEnergy = 0.5;

    arma::rowvec runVector;
    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.1, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.2, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.3, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.4, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.6, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.7, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.8, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.9, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);
}


TEST_CASE("Local energy, 10 particles, 1D, alpha != 0.5", "[Hamiltonian, alpha not equal to 0.5, 10 particles, 1 dimension]")
{
    int nParticles = 10;
    int nDimensions = 1;
    int nCycles = 1e6;
    double stepLength = 0.1;
    double dt = 0.01;
    int cycleStepToFile = 1e6;

    double energy = 0;
    double localEnergy = 5;

    arma::rowvec runVector;
    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.1, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.2, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.3, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.4, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.6, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.7, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.8, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.9, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);
}


TEST_CASE("Local energy, 1 particle, 2D, alpha != 0.5", "[Hamiltonian, alpha not equal to 0.5, 1 particle, 2 dimensions]")
{
    int nParticles = 1;
    int nDimensions = 2;
    int nCycles = 1e6;
    double stepLength = 0.1;
    double dt = 0.01;
    int cycleStepToFile = 1e6;

    double energy = 0;
    double localEnergy = 1;

    arma::rowvec runVector;
    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.1, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.2, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.3, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.4, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.6, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.7, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.8, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.9, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);
}


TEST_CASE("Local energy, 10 particles, 2D, alpha != 0.5", "[Hamiltonian, alpha not equal to 0.5, 10 particles, 2 dimensions]")
{
    int nParticles = 10;
    int nDimensions = 2;
    int nCycles = 1e6;
    double stepLength = 0.1;
    double dt = 0.01;
    int cycleStepToFile = 1e6;

    double energy = 0;
    double localEnergy = 10;

    arma::rowvec runVector;
    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.1, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.2, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.3, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.4, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.6, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.7, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.8, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.9, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);
}


TEST_CASE("Local energy, 1 particle, 3D, alpha != 0.5", "[Hamiltonian, alpha not equal to 0.5, 1 particle, 3 dimensions]")
{
    int nParticles = 1;
    int nDimensions = 3;
    int nCycles = 1e6;
    double stepLength = 0.1;
    double dt = 0.01;
    int cycleStepToFile = 1e6;

    double energy = 0;
    double localEnergy = 1.5;

    arma::rowvec runVector;
    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.1, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.2, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.3, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.4, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.6, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.7, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.8, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.9, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);
}


TEST_CASE("Local energy, 10 particles, 3D, alpha != 0.5", "[Hamiltonian, alpha not equal to 0.5, 10 particles, 3 dimensions]")
{
    int nParticles = 10;
    int nDimensions = 3;
    int nCycles = 1e6;
    double stepLength = 0.1;
    double dt = 0.01;
    int cycleStepToFile = 1e6;

    double energy = 0;
    double localEnergy = 15;

    arma::rowvec runVector;
    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.1, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.2, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.3, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.4, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.6, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.7, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.8, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);

    runVector = VMC->RunMonteCarloIntegration(nParticles, nDimensions, nCycles, 0.9, stepLength, dt, cycleStepToFile);
    energy = runVector(1);
    REQUIRE(energy > localEnergy);
}

#endif // CATCH_HPP
