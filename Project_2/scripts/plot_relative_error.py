import numpy as np
import matplotlib.pyplot as plt

analytical_data = "../results/energies_analytic.txt"
numerical_data = "../results/energies_numerical.txt"


def RelativeError(data):
    cycles, analytical_energies = np.loadtxt(data, unpack=True)
    exact_energy = 2.0
    absolute_error = abs(analytical_energies - exact_energy)
    relative_error = absolute_error / exact_energy
    x = cycles[:2000]
    y = relative_error[:2000]
    return x, y


xa, ya = RelativeError(analytical_data)
xn, yn = RelativeError(numerical_data)

plt.plot(xa, ya, color="red", linestyle="-", label="Analytical")
plt.plot(xn, yn, color="blue", linestyle="-.", label="Numerical")
plt.xlabel("Cycles")
plt.ylabel(r"Energy $[\hbar \omega]$")
plt.legend(loc="best")
plt.savefig("../results/rel_error.png")
plt.show()
