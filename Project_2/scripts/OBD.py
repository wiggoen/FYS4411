# Python 3.6
# One-body density

import matplotlib.pyplot as plt
import numpy as np

Int_J_path = "../results/histogram-Int-J-2p-2up27.txt"
Int_noJ_path = "../results/histogram-Int-noJ-2p-2up27.txt"
noInt_J_path = "../results/histogram-noInt-J-2p-2up27.txt"
noInt_noJ_path = "../results/histogram-noInt-noJ-2p-2up27.txt"

Int_J = np.loadtxt(Int_J_path)
Int_noJ = np.loadtxt(Int_noJ_path)
noInt_J = np.loadtxt(noInt_J_path)
noInt_noJ = np.loadtxt(noInt_noJ_path)

distance = np.linspace(0, 4, len(Int_J)+1)
norm = distance[1:]**2/(2*np.pi*distance[1:]*distance[1] - np.pi*distance[1]*distance[1])

Int_J = Int_J*norm
Int_noJ = Int_noJ*norm
noInt_J = noInt_J*norm
noInt_noJ = noInt_noJ*norm

fig = plt.figure()
plt.plot(distance[1:], Int_J/sum(Int_J), linestyle='None', marker='o',
         markersize=0.6, label="Interaction on, Jastrow on")
plt.plot(distance[1:], Int_noJ/sum(Int_noJ), linestyle='None', marker='o',
         markersize=0.6, label="Interaction on, Jastrow off")
plt.plot(distance[1:], noInt_J/sum(noInt_J), linestyle='None', marker='o',
         markersize=0.6, label="Interaction off, Jastrow on")
plt.plot(distance[1:], noInt_noJ/sum(noInt_noJ), linestyle='None', marker='o',
         markersize=0.6, label="Interaction off, Jastrow off")

plt.xlabel(r"$|\mathbf{r}|$", fontsize=15)
plt.ylabel(r"$\rho(\mathbf{r})$", fontsize=15)
plt.legend(loc=0, fontsize=15)
fig.set_tight_layout(True)  # Minimizing overlap of labels
plt.savefig("../plots/obd_2p.png")
plt.show()
