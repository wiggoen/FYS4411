# Python 3.6
# One-body density

import matplotlib.pyplot as plt
import numpy as np

int_J_data = "../results/histogram-int-J-2p.txt"
int_noJ_data = "../results/histogram-int-noJ-2p.txt"
noInt_J_data = "../results/histogram-non-int-J-2p.txt"
noInt_noJ_data = "../results/histogram-non-int-noJ-2p.txt"

int_J = np.loadtxt(int_J_data)
int_noJ = np.loadtxt(int_noJ_data)
noInt_J = np.loadtxt(noInt_J_data)
noInt_noJ = np.loadtxt(noInt_noJ_data)

distance_int_J = np.linspace(0, 2, len(int_J))
distance_int_noJ = np.linspace(0, 2, len(int_noJ))
distance_noInt_J = np.linspace(0, 2, len(noInt_J))
distance_noInt_noJ = np.linspace(0, 2, len(noInt_noJ))


fig = plt.figure()
plt.plot(distance_int_J, int_J/sum(int_J),
         linestyle='None', marker='o', markersize=0.6,
         label="Interaction on, Jastrow on")
plt.plot(distance_int_noJ, int_noJ/sum(int_noJ),
         linestyle='None', marker='o', markersize=0.6,
         label="Interaction on, Jastrow off")
plt.plot(distance_noInt_J, noInt_J/sum(noInt_J),
         linestyle='None', marker='o', markersize=0.6,
         label="Interaction off, Jastrow on")
plt.plot(distance_noInt_noJ, noInt_noJ/sum(noInt_noJ),
         linestyle='None', marker='o', markersize=0.6,
         label="Interaction off, Jastrow off")
plt.xlabel(r"$|\mathbf{r}|$", fontsize=15)
plt.ylabel(r"$\rho(\mathbf{r})$", fontsize=15)
plt.legend(loc=0, fontsize=15)
fig.set_tight_layout(True)  # Minimizing overlap of labels
plt.savefig("../plots/obd_2p.png")
plt.show()
