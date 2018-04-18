# Script for Python 3

import subprocess
import os


def run_cmd(cmd):
    print(" ")
    print(">> {}".format(cmd))
    subprocess.call(cmd, shell=True)


# File names
fname1 = "../main"
fname2 = "../src/hamiltonian"
fname3 = "../src/tests"
fname4 = "../src/variationalmontecarlo"
fname5 = "../src/wavefunction"


# Parameters
nParticles = 6
nDimensions = 3
nCycles = 1e6
alpha = 0.5
stepLength = 0.1
timeStep = 0.01
cycleStepToFile = nCycles
trials = 1
# Choices: "bf" (Brute force) or "im" (Importance)
samplingType = "im"
# Choices: "ana" (Analytical) or "num" (Numerical) or "int" (Interaction)
derivationType = "int"


# Library used
library = "-L/usr/local/Cellar/armadillo/8.400.0/lib/"

# Location of include files / header files
includes = "-I/usr/local/Cellar/armadillo/8.400.0/include/ -I/Users/trondwj/GitHub/FYS4411/Project_1"


# Compile
run_cmd("g++ -std=c++11 -O3 {} {} -o {}.out {}.cpp {}.cpp {}.cpp {}.cpp {}.cpp -larmadillo".format(library, includes, fname1, fname1, fname2, fname3, fname4, fname5))

# Run
run_cmd("./{}.out {} {} {} {} {} {} {} {} {} {}".format(fname1, nParticles, nDimensions, nCycles, alpha, stepLength, timeStep, cycleStepToFile, trials, samplingType, derivationType))
