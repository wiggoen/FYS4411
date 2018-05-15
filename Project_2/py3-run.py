# Script for Python 3
import os

machine = "Linux"
#machine = "Mac"


def run_cmd(cmd):
    print(" ")
    print(">> {}".format(cmd))
    os.system(cmd)


# Compile path
src = "src/*.cpp"

# Parameter file
parameterfile = "parameters.json"

# Library used
if (machine == "Linux"):
    library = "-L/usr/lib"
elif (machine == "Mac"):
    library = "-L/usr/local/Cellar/armadillo/8.400.0/lib/"

# Location of include files / header files
if (machine == "Linux"):
    includes = "-I/usr/include -I/home/twj/Documents/GitHub/FYS4411/Project_2"
elif (machine == "Mac"):
    includes = "-I/usr/local/Cellar/armadillo/8.400.0/include/ -I/Users/trondwj/GitHub/FYS4411/Project_2"

# Compile
run_cmd("g++ -std=c++11 -O3 {} {} -o main.out main.cpp {} -larmadillo".format(library, includes, src))

# Run
run_cmd("./main.out {}".format(parameterfile))
