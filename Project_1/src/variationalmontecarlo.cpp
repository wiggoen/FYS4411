#include "inc/variationalmontecarlo.h"
#include <random>
#include <iomanip>
#include <iostream>

variationalMonteCarlo::variationalMonteCarlo()
{
    // Propose a new position R by moving one boson at the time
    // Calculate new psi
    // Pick random number r in [0,1]
    // Test if r is smaller or equal to |psi_T(R')|^2/|psi_T(R')|^2 ??
    // If yes: accept new position
        // Calculate delta E_L(R)
        // Update E = E + delta E_L
        // E^2 = E^2 + (delta E_L)^2
}

variationalMonteCarlo::~variationalMonteCarlo()
{

}

