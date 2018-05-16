# Script for Python 3
import os
import json


# Choice of machine type
#machine = "Linux"
machine = "Mac"


def run_cmd(cmd):
    print(" ")
    print(">> {}".format(cmd))
    os.system(cmd)


# Compile path
src = "src/*.cpp"

# Parameter file
parameterfile = "parameters.json"
with open(parameterfile) as data:
    parameters = json.load(data)
UseMPI = parameters["Use MPI"]
print(UseMPI)

# Library used
if (machine == "Linux"):
    library = "-L/usr/lib"
elif (machine == "Mac"):
    library = "-L/usr/local/Cellar/armadillo/8.500.1/lib/ -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/3.1.0/lib -lmpi"

# Location of include files / header files
if (machine == "Linux"):
    includes = "-I/usr/include -I/home/twj/Documents/GitHub/FYS4411/Project_2"
elif (machine == "Mac"):
    includes = "-I/usr/local/Cellar/armadillo/8.500.1/include/ -I/Users/trondwj/GitHub/FYS4411/Project_2 -I/usr/local/Cellar/open-mpi/3.1.0/include"

# Compile
if (UseMPI):
    run_cmd("mpicxx -std=c++11 -O3 {} {} -o main.out main.cpp {} -larmadillo".format(library, includes, src))
else:
    run_cmd("g++ -std=c++11 -O3 {} {} -o main.out main.cpp {} -larmadillo".format(library, includes, src))

# Run
if (UseMPI):
    run_cmd("mpirun -n 2 ./main.out {}".format(parameterfile))
else:
    run_cmd("./main.out {}".format(parameterfile))
