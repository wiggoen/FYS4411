#define CATCH_CONFIG_RUNNER // Configure Catch to use this main, and not its own.
#include "inc/catch.hpp"
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
        return 0;
    }
}
