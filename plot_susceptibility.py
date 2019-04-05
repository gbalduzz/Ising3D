#!/usr/bin/python3

from matplotlib import pyplot as plt
import numpy as np
import json

file = json.load(open("input.json"))
Ls = file["Ls"]

def jackKnife(procedure, x, n):
    block_size = int(np.ceil(len(x) / n))
    f = []
    for i in range(n):
        start = int(i*block_size)
        end = int(min((i+1)*block_size, len(x)-1))
        x_loo = np.concatenate([x[:start], x[end:]])
        f.append(procedure(x_loo))

    f_bar = procedure(x)
    f = np.array(f)

    delta = np.sqrt((n-1)/n * np.sum((f-f_bar)**2))
    estimate = n * f_bar - (n-1) * np.mean(f)
    return estimate, delta

def sort(array, col):
    return array[array[:,col].argsort()]

for L in Ls:
    data = np.loadtxt("outputs/magnetization_L"+str(L)+".txt")
    betas = data[:, 0]
    out = []

    for id, beta in enumerate(betas):

	# Note: the c code outputs a signed magnetization density, here we need the total, hence we
	#       multiply by the volume.
        M = np.abs(data[id, 1:]) * L**3

        def susceptibility(m) :
            m2 = (m**2).mean()
            m = m.mean()
            return beta*(m2-m**2)

        val, err = jackKnife(susceptibility, M, 20)
        out.append([1./beta, val, err])

    out = sort(np.array(out), 0)
    plt.errorbar(out[:,0], out[:,1], fmt = "--o", yerr = out[:,2], label = "L="+str(L))

plt.legend(loc="best")
plt.xlabel(r"$T$")
plt.ylabel(r"$\chi$")
plt.title("Magnetic susceptibility.")
plt.show()
