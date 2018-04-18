# Script for Python 2

import subprocess
import os


def run_cmd(cmd):
    print ">> ", cmd
    subprocess.call(cmd, shell=True)


# File names
fname1 = "../main"
fname2 = "../src/hamiltonian"
fname3 = "../src/tests"
fname4 = "../src/variationalmontecarlo"
fname5 = "../src/wavefunction"


# Parameters
nParticles = [1, 10, 100, 500]
nDimensions = [1, 2, 3]
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
library = "-L/usr/include"

# Location of include files / header files
includes = "-I/usr/local/include -I/home/line/github/FYS4411/Project_1"


# Compile
run_cmd("g++ -std=c++11 -O3 %s %s -o %s.out %s.cpp %s.cpp %s.cpp %s.cpp %s.cpp -larmadillo" % (library, includes, fname1, fname1, fname2, fname3, fname4, fname5))

# Run
for i in range(len(nDimensions)):
    for j in range(len(nParticles)):
        run_cmd("./%s.out %s %s %s %s %s %s %s %s %s %s" % (fname1, nParticles[j], nDimensions[i], nCycles, alpha, stepLength, timeStep, cycleStepToFile, trials, samplingType, derivationType))
