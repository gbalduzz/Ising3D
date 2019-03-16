#!/usr/bin/python3

from matplotlib import pyplot as plt
from sys import argv
import numpy as np
import glob
import re

if len(argv) < 2:
	print("Usage: pythhon binder_cumulant.py <data folder>")
	exit(-1)

folder = argv[1]


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

def binderCumulant(x) :
    x2 = (x**2).mean()
    x4 = (x**4).mean()
    return 1 - x4/(3*x2**2)

def sort(array, col):
    return array[array[:,col].argsort()]


# Find the length of the simulated clusters.
files = glob.glob(folder + "/magnetization_L*" + ".txt")
Ls = []
for file in files:
	L = re.search("[0-9]+", file).group(0)
	Ls.append(L)

for L in Ls:
    data = np.loadtxt(folder + "/magnetization_L"+str(L)+".txt")
    betas = data[:, 0]
    out = []

    for id, beta in enumerate(betas):

        M = data[id, 1:]
        b, db = jackKnife(binderCumulant, M, 20)

        out.append([1./beta, b, db])

    out = sort(np.array(out), 0)
    np.savetxt(folder + "/binder_cumulants_L" + str(L) + ".txt", out, header='T\tU\tdelta U')
    plt.errorbar(out[:,0], out[:,1], fmt = "--o", yerr = out[:,2], label = "L="+str(L))

plt.legend(loc="best")
plt.xlabel(r"$T$")
plt.ylabel(r"$U$")
plt.title("Binder cumulant.")

plt.savefig("binder_cumulant.pdf", format = 'pdf')
plt.show()
