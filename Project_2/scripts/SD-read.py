import numpy as np

data = "../results/SD-run.txt"
SD_data = np.loadtxt(data, skiprows=1)

minima = np.amin(SD_data[:, 3])
print("Energy minima = {}".format(minima))
