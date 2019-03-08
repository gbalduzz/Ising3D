#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt("output.txt")

beta = data[:,0]
E = data[:,1]
M = data[:,2]
E2 = data[:,3]
M2 = data[:,4]

c = beta**2 * (E2 - E**2)
chi = beta * (M2 - M**2)

plt.figure(0)
plt.plot(beta, chi)
plt.xlabel("1/T")
plt.ylabel(r"$\chi$")

plt.figure(1)
plt.plot(beta, c)
plt.xlabel("1/T")
plt.ylabel(r"$c$")

plt.figure(2)
energies = np.loadtxt("energies.txt")
id = np.linspace(1, energies.shape[0], energies.shape[0])
plt.plot(id, energies)
plt.xlabel("sweep")
plt.ylabel("E")

plt.show()
