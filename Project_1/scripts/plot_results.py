import numpy as np
import matplotlib.pyplot as plt

# Choose 2 for Python 2.7
# or 3 for Python 3.x
Python_version = 3

if Python_version == 2:
	# Python 2.7
	cyclestep = []
	energy    = []

	f = open("results.txt","r")
	lines = f.readlines()

	for line in lines:
		cyclestep.append(line.split('     ')[0])
		energy.append(line.split('     ')[1])
	f.close()

elif Python_version == 3:
	# Python 3.6
	datafile = "results.txt"
	cyclestep = np.loadtxt(datafile, usecols=0, unpack=True)
	energy = np.loadtxt(datafile, usecols=1, unpack=True)


fig = plt.figure()
plt.plot(cyclestep,energy,linewidth=1)
plt.xlabel("Cycles",fontsize=15)
plt.ylabel("Energy",fontsize=15)
#plt.legend(fonsize=15)
#plt.tick_params(axis='x', labelsize=12)
#plt.tick_params(axis='y', labelsize=12)
fig.set_tight_layout(True)  # Minimizing overlap of labels
plt.show()
