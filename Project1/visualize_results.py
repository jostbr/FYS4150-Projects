
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

problem_sizes = np.logspace(1, 4, 4, endpoint = True).astype(int)
problem_sizes = np.arange(10, 101, 10)
num_boundary_points = 2

plt.figure(figsize = (8, 7))
plt.style.use("seaborn")

for size in problem_sizes:
    full_system_size = size + num_boundary_points     # Account for boundary points
    x = np.zeros(full_system_size)
    numerical = np.zeros(full_system_size)
    analytical = np.zeros(full_system_size)
    #rel_error = np.zeros(size)

    with open("result{}.txt".format(size), "r") as results:
        for line_num, line_string in enumerate(results):
            if (line_num > 0):
                words = line_string.split()
                x[line_num-1] = float(words[0])
                numerical[line_num-1] = float(words[1])
                analytical[line_num-1] = float(words[2])
                #rel_error[line_num-1] = float(words[3])

        plt.plot(x, numerical, linewidth = 2, label = "Numerical, N = {}".format(size))

plt.plot(x, analytical, linewidth = 2, label = "Analytical")
plt.xlabel("$x$-coordinate", fontname = "serif", fontsize = 16)
plt.ylabel("$u(x)$", fontname = "serif", fontsize = 16)
plt.title("Solution to 1D Poisson eq. (gen. alg.)",
    fontname = "serif", fontsize = 18)
plt.legend()
plt.show()
