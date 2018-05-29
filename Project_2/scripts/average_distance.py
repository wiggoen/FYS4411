import numpy as np

data = "../results/distances.txt"
distance = np.loadtxt(data, usecols=1)

average_distance = np.mean(distance)
print("Average distance = {}".format(average_distance))
