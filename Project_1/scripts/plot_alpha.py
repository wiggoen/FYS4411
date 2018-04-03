import numpy as np
import matplotlib.pyplot as plt

filename = open("1dim.txt", "r")
lines = filename.readlines()

energy = np.loadtxt("1dim.txt")
print energy 

final_energies = []


sims = 0
counter = 0
for i in range(len(energy)):
	sims += energy[i]
	counter += 1
	print i	
	if (counter == 10): 
		final_energies.append(sims/10.0)
		sims = 0
		counter = 0

print final_energies 

print len(final_energies)

alpha = np.linspace(0.1,0.9,9)

#print alpha
#print final_energies[0:9:1]
#print len(final_energies[10:19:1])
#print len(final_energies[20:29:1])
#print len(final_energies[10:19:1])
#plt.plot(alpha,final_energies[20:29:1])
fig = plt.figure()

plt.plot(alpha,np.log(final_energies[0:9:1]),label="1D, 1 particle")
plt.plot(alpha,np.log(final_energies[9:18:1]),label="1D, 10 particles")
plt.plot(alpha,np.log(final_energies[18:27:1]),label="1D, 100 particles")
plt.ylabel("Log(Energy) [Relative units]",fontsize=15)
plt.legend(fontsize=15)
plt.xlabel("Alpha",fontsize=15)
plt.tick_params(axis='x', labelsize=12)
plt.tick_params(axis='y', labelsize=12)
fig.set_tight_layout(True)

plt.show()
