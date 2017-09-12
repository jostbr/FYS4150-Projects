
# Script for plotting the numerical vs. the analytical solution.
# Depends on having generated result10,txt, result100.txt and result1000.txt.

# You might have to change file path for result files depending on your build settings.

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

problem_sizes = np.logspace(1, 3, 3, endpoint = True).astype(int)
#problem_sizes = np.arange(10, 101, 10)
#problem_sizes = np.arange(7, 84, 7)
num_boundary_points = 2

plt.figure(figsize = (8, 7))
plt.style.use("seaborn")

for size in problem_sizes:
    full_system_size = size + num_boundary_points     # Account for boundary points
    x = np.zeros(full_system_size)
    numerical = np.zeros(full_system_size)
    analytical = np.zeros(full_system_size)

    with open("../build-Project1/result{}.txt".format(size), "r") as results:
        for line_num, line_string in enumerate(results):
            if (line_num > 0):
                words = line_string.split()
                x[line_num-1] = float(words[0])
                numerical[line_num-1] = float(words[1])
                analytical[line_num-1] = float(words[2])

        plt.plot(x, numerical, linewidth = 2,
            label = "Numerical: N = {}, h = {:.2e}".format(full_system_size, 1/(full_system_size-1)))

plt.plot(x, analytical, linewidth = 2, label = "Analytical: $u(x)=1-(1-e^{-10})x-e^{-10x}$")
plt.xlabel("$x$-coordinate [dim-less]", fontname = "serif", fontsize = 15)
plt.ylabel("$u(x)$ [dim-less]", fontname = "serif", fontsize = 15)
plt.title("Solution to 1D Poisson equation, source: $g(x)=100e^{-10x}$",
    fontname = "serif", fontsize = 16)
plt.legend()
plt.show()
