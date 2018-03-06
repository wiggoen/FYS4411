#ifndef CATCH_HPP
#include "inc/catch.hpp"
#include "inc/variationalmontecarlo.h"


VariationalMonteCarlo *VMC = new VariationalMonteCarlo();

TEST_CASE("Local energy", "[Hamiltonian]")
{
    double energy = 0;

    energy = VMC->RunMonteCarloIntegration(1, 1, 1e6, 0.5, 0.1, 1e6);
    REQUIRE(energy == 0.5);

    energy = VMC->RunMonteCarloIntegration(1, 2, 1e6, 0.5, 0.1, 1e6);
    REQUIRE(energy == 1.0);

    energy = VMC->RunMonteCarloIntegration(1, 3, 1e6, 0.5, 0.1, 1e6);
    REQUIRE(energy == 1.5);

    energy = VMC->RunMonteCarloIntegration(10, 1, 1e6, 0.5, 0.1, 1e6);
    REQUIRE(energy == 5);

    energy = VMC->RunMonteCarloIntegration(10, 2, 1e6, 0.5, 0.1, 1e6);
    REQUIRE(energy == 10);

    energy = VMC->RunMonteCarloIntegration(10, 3, 1e6, 0.5, 0.1, 1e6);
    REQUIRE(energy == 15);
}

#endif // CATCH_HPP
