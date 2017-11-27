
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def get_data(filename):
    data = np.loadtxt(filename, dtype = np.float64)

    t = data[:, 0]
    psi = data[:, 1:]
    x = np.linspace(0, 1, 41)
    y = np.linspace(0, 1, 41)

    return x, y, t, psi

def contour_results(x, y, psi):
	fig = plt.figure(figsize = (9, 7))
	plt.title("Contour field of $\psi$", fontname = "serif", fontsize = 17)
	plt.xlabel("x [dim-less]", fontname = "serif", fontsize = 12)
	plt.ylabel("y [dim-less]", fontname = "serif", fontsize = 12)
	CS = plt.contourf(x, y, psi.transpose(), 20, cmap = plt.cm.RdBu_r)
	plt.colorbar(CS, orientation = "vertical")


filename_2D = "../build-Project5/results_2d.txt"
x, y, t, psi = get_data(filename_2D)

psi_0 = psi[0, :].reshape((41, 41))

contour_results(x, y, psi_0)

plt.show()