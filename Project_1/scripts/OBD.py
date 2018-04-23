# Python 3.6
# One-body density

import matplotlib.pyplot as plt
import numpy as np

NoInteractionData = "../results/histogram-non-interacting-20p-alpha0495.txt"
InteractionData = "../results/histogram-interacting-20p-alpha0495.txt"
NoInteraction = np.loadtxt(NoInteractionData)
Interaction = np.loadtxt(InteractionData)

distance_NoInteraction = np.linspace(0, 2, len(NoInteraction))
distance_Interaction = np.linspace(0, 2, len(Interaction))

fig = plt.figure()
plt.plot(distance_NoInteraction, NoInteraction/sum(NoInteraction),
         linestyle='None', marker='o', markersize=0.6,
         label="No Interaction")
plt.plot(distance_Interaction, Interaction/sum(Interaction),
         linestyle='None', marker='o', markersize=0.6,
         label="Interaction")
plt.xlabel(r"$|\mathbf{r}|$", fontsize=15)
plt.ylabel(r"$\rho(\mathbf{r})$", fontsize=15)
plt.legend(loc=0, fontsize=15)
fig.set_tight_layout(True)  # Minimizing overlap of labels
plt.savefig("../plots/obd_20p.png")
plt.show()
