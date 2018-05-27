# Script for Python 3
import os
import sys
import json
import subprocess

c_compile_sh = subprocess.run("mpicc --showme:compile", stdout=subprocess.PIPE, shell=True)
c_compile = c_compile_sh.stdout.decode("utf-8").strip()
cxx_link_sh = subprocess.run("mpicxx --showme:link", stdout=subprocess.PIPE, shell=True)
cxx_link = cxx_link_sh.stdout.decode("utf-8").strip()
cxx_compile_sh = subprocess.run("mpicxx --showme:compile", stdout=subprocess.PIPE, shell=True)
cxx_compile = cxx_compile_sh.stdout.decode("utf-8").strip()


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
#print(UseMPI)

# Location of include files / header files and libraries used
if sys.platform.startswith("linux"):    # Linux specific
    armadillo_include = "-I/usr/include"
    project_folder = "-I/home/twj/Documents/GitHub/FYS4411/Project_2"
    if UseMPI:
        includes = "{} {} -DMPICH_IGNORE_CXX_SEEK {} {}".format(c_compile, cxx_compile, armadillo_include, project_folder)
        libraries = "{} -L/usr/lib".format(cxx_link)
    else:
        includes = "{} {}".format(armadillo_include, project_folder)
        libraries = "-L/usr/lib"
elif sys.platform.startswith("darwin"): # Mac specific
    armadillo_include = "-I/usr/local/Cellar/armadillo/8.500.1/include/"
    project_folder = "-I/Users/trondwj/GitHub/FYS4411/Project_2"
    if UseMPI:
        includes = "{} {} -DMPICH_IGNORE_CXX_SEEK {} {}".format(c_compile, cxx_compile, armadillo_include, project_folder)
        libraries = "{} -L/usr/local/Cellar/armadillo/8.500.1/lib/".format(cxx_link)
    else:
        includes = "{} {}".format(armadillo_include, project_folder)
        libraries = "-L/usr/local/Cellar/armadillo/8.500.1/lib/"

# Compile
if (UseMPI):
    run_cmd("mpicxx -std=c++11 -O3 {} {} -o main.out main.cpp {} -larmadillo -DMPI_ON".format(includes, libraries, src))
else:
    run_cmd("g++ -Wall -std=c++11 -O3 {} {} -o main.out main.cpp {} -llapack -lblas -larmadillo".format(includes, libraries, src))

# Run
if (UseMPI):
    run_cmd("mpirun -n 2 ./main.out {}".format(parameterfile))
else:
    run_cmd("./main.out {}".format(parameterfile))
