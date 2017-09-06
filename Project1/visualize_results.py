
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

problem_sizes = np.logspace(1, 3, 3, endpoint = True).astype(int)
#problem_sizes = np.arange(10, 101, 10)
num_boundary_points = 2

plt.figure(figsize = (8, 7))
plt.style.use("ggplot")

for size in problem_sizes:
    full_system_size = size + num_boundary_points     # Account for boundary points
    x = np.zeros(full_system_size)
    numerical = np.zeros(full_system_size)
    analytical = np.zeros(full_system_size)
    #rel_error = np.zeros(size)

    with open("../build-Project1/result{}.txt".format(size), "r") as results:
        for line_num, line_string in enumerate(results):
            if (line_num > 0):
                words = line_string.split()
                x[line_num-1] = float(words[0])
                numerical[line_num-1] = float(words[1])
                analytical[line_num-1] = float(words[2])
                #rel_error[line_num-1] = float(words[3])

        plt.plot(x, numerical, linewidth = 2,
            label = "Numerical, N = {}, h = {:.2e}".format(full_system_size, 1/(full_system_size-1)))

plt.plot(x, analytical, linewidth = 2, label = "Analytical")
plt.xlabel("$x$-coordinate [dim-less]", fontname = "serif", fontsize = 16)
plt.ylabel("$u(x)$ [dim-less]", fontname = "serif", fontsize = 16)
plt.title("Solution to 1D Poisson equation",
    fontname = "serif", fontsize = 18)
plt.legend()
plt.show()
