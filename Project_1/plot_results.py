import numpy as np
import matplotlib.pyplot as plt


cyclestep = []
energy    = []

f = open("results.txt","r")
lines = f.readlines()

for line in lines:
	cyclestep.append(line.split('    ')[0])
	energy.append(line.split('    ')[1])
f.close()

plt.plot(cyclestep,energy,linewidth=1)
plt.xlabel("Cycles",fontsize=15)
plt.ylabel("Energy",fontsize=15)
#plt.legend(fonsize=15)
#plt.tick_params(axis='x', labelsize=12)
#plt.tick_params(axis='y', labelsize=12)
plt.show()
