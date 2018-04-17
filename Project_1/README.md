# Variational Monte Carlo of Bosonic systems

In collaboration with [Line](https://github.com/linegpe).

## Ways of running the program (you can change the parameters given by default):
- Python scripts (you have to change to the library/include paths to fit your installation of [Armadillo](http://arma.sourceforge.net))
  - [Python 2](scripts/run.py) [relative link]
  - [Python 3](scripts/py3-run.py) [relative link]
- Manually by running it in [Qt Creator](https://www.qt.io) (you have to link [Armadillo](http://arma.sourceforge.net) to the .pro-file)
- Manually by using the command line with command line arguments (basically the same as the Python scripts provide)


## Running automatic tests of the code
This is unfortunately done manually. In 'tests.cpp' you can choose what to test:
```
// CHOOSE SAMPLING METHOD                       <<< --- CHOOSE ONLY ONE FOR TESTING
std::string samplingType = "BruteForce";
//std::string samplingType = "Importance";

// CHOOSE INTEGRATION METHOD                    <<< --- CHOOSE ONLY ONE FOR TESTING
std::string derivationType = "Analytical";
//std::string derivationType = "Numerical";
//std::string derivationType = "Interaction";
```
and in top of 'main.cpp' you have to change
```
#define TEST false
```
into
```
#define TEST true
```
Remember to put it back to 'false' when running the program!