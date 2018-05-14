# Script for Python 3
import os


def run_cmd(cmd):
    print(" ")
    print(">> {}".format(cmd))
    os.system(cmd)


# Compile path
src = "src/*.cpp"

# Parameter file
parameterfile = "parameters.json"

# Library used
library = "-L/usr/local/Cellar/armadillo/8.400.0/lib/"

# Location of include files / header files
includes = "-I/usr/local/Cellar/armadillo/8.400.0/include/ -I/Users/trondwj/GitHub/FYS4411/Project_2"

# Compile
run_cmd("g++ -std=c++11 -O3 {} {} -o main.out main.cpp {} -larmadillo".format(library, includes, src))

# Run
run_cmd("./main.out {}".format(parameterfile))
