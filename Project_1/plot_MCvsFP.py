import numpy as np
import matplotlib.pyplot as plt

MC5 = np.loadtxt("MCa05p100d3.txt", unpack=True)
MC4 = np.loadtxt("MCa04p100d3.txt", unpack=True)
MC6 = np.loadtxt("MCa06p100d3.txt", unpack=True)

FP5 = np.loadtxt("FPa05p100d3.txt", unpack=True)
FP4 = np.loadtxt("FPa04p100d3.txt", unpack=True)
FP6 = np.loadtxt("FPa06p100d3.txt", unpack=True)







plt.plot(MC5[0],MC5[1],label = "MC,a=0.5")
plt.plot(MC4[0],MC4[1],label = "MC,a=0.4")
plt.plot(MC6[0],MC6[1],label = "MC,a=0.6")
plt.plot(FP5[0],FP5[1],label = "FP,a=0.5")
plt.plot(FP4[0],FP4[1],label = "FP,a=0.4")
plt.plot(FP6[0],FP6[1],label = "FP,a=0.6")
plt.legend(fontsize=15)
plt.axis([-1000,1000000,145,160])
plt.ylabel("Energy", fontsize=15)
plt.xlabel("Cycles", fontsize=15)
plt.show()

