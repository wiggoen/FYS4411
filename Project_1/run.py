import subprocess
import os

def run_cmd(cmd):
    print ">> ", cmd
    subprocess.call(cmd,shell=True)

fname1 = "main"
fname2 = "src/hamiltonian"
fname3 = "src/tests"
fname4 = "src/variationalmontecarlo"
fname5 = "src/wavefunction"

nParticles = [1,10,100,500]
nDim = [1,2,3]
nCycles = 1e6
alpha = 0.5
stepLength = 0.1
nDataPoints = 1e6

run_cmd("g++ -std=c++11 -o2 -L/usr/include -I/usr/local/include -I/home/line/github/FYS4411/Project_1 -o %s.out %s.cpp %s.cpp %s.cpp %s.cpp %s.cpp -larmadillo" % (fname1, fname1, fname2, fname3, fname4, fname5))


for i in range(len(nDim)):
    for j in range(len(nParticles)):
		run_cmd("./%s.out %s %s %s %s %s %s" % (fname1, nParticles[j], nDim[i], nCycles, alpha, stepLength, nDataPoints))
