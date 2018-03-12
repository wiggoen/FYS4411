import numpy as np
import matplotlib.pyplot as plt

a1p1d = np.loadtxt("1p1d.txt", unpack=True, skiprows=1)
a1p2d = np.loadtxt("1p2d.txt", unpack=True, skiprows=1)
a1p3d = np.loadtxt("1p3d.txt", unpack=True, skiprows=1)
a10p1d = np.loadtxt("10p1d.txt", unpack=True, skiprows=1)
a10p2d = np.loadtxt("10p2d.txt", unpack=True, skiprows=1)
a10p3d = np.loadtxt("10p3d.txt", unpack=True, skiprows=1)
a100p1d = np.loadtxt("100p1d.txt", unpack=True, skiprows=1)
a100p2d = np.loadtxt("100p2d.txt", unpack=True, skiprows=1)
a100p3d = np.loadtxt("100p3d.txt", unpack=True, skiprows=1)
a500p1d = np.loadtxt("500p1d.txt", unpack=True, skiprows=1)
a500p2d = np.loadtxt("500p2d.txt", unpack=True, skiprows=1)
#a500p3d = np.loadtxt("500p3d.txt", unpack=True, skiprows=1)


plt.plot(a1p1d[3],np.log(a1p1d[6]),label="1p,1D")
plt.plot(a1p1d[3],np.log(a1p2d[6]),label="1p,2D")
plt.plot(a1p1d[3],np.log(a1p3d[6]),label="1p,3D")
plt.plot(a1p1d[3],np.log(a10p1d[6]),label="10p,1D")
plt.plot(a1p1d[3],np.log(a10p2d[6]),label="10p,2D")
plt.plot(a1p1d[3],np.log(a10p3d[6]),label="10p,3D")
plt.plot(a1p1d[3],np.log(a100p1d[6]),label="100p,1D")
plt.plot(a1p1d[3],np.log(a100p2d[6]),label="100p,2D")
plt.plot(a1p1d[3],np.log(a100p3d[6]),label="100p,3D")
plt.plot(a1p1d[3],np.log(a500p1d[6]),label="500p,1D")
plt.plot(a1p1d[3],np.log(a500p2d[6]),label="500p,2D")
#plt.plot(a1p1d[3],np.log(a500p3d[6]),label="500p,2D")
plt.xlabel("Alpha",fontsize=15)
plt.ylabel("Log(energy) [Relative units]",fontsize=15)
plt.legend(fontsize=15)
plt.show()

