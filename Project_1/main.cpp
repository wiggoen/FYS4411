#define CATCH_CONFIG_RUNNER // Configure Catch to use this main, and not its own.
#include "inc/catch.hpp"
#include "inc/matrix.h"
#include <iostream>

#define TEST false // Change to true when testing and to false when running the program.

int runCatchTests() {
    return Catch::Session().run();
}

int main() {
    if (TEST) {
        std::cout << "Running tests" << std::endl;
        return runCatchTests();
    } else {

        std::cout << "Hello World!" << std::endl;
        Matrix* sys = new Matrix(5, 3);

        double beta = 1;
        int N = 1;


        // Timing the algorithm
        clock_t start, finish;
        start = clock();

        // Timing finished
        finish = clock();
        double time = ( double (finish - start)/CLOCKS_PER_SEC);

        std::cout << std::endl << "Run time: " << time << " sec.";

        return 0;
    }


}
