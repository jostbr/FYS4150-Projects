
# Script for plotting the numerical Armadillo solution vs. the analytical.
# Depends on having generated results with the armadillo solve() algorithm.

# You might have to change file path for result files depending on your build settings.

import numpy as np
import matplotlib.pyplot as plt

x = list()
numerical = list()
analytical = list()

with open("../build-Project1/arma_results.txt", "r") as result_file:
    for line_num, line_string in enumerate(result_file):
        if (line_num > 0):
            words = line_string.split()
            x.append(float(words[0]))
            numerical.append(float(words[1]))
            analytical.append(float(words[2]))

plt.figure(figsize = (8, 7))
plt.style.use("seaborn")
plt.plot(x, numerical, linewidth = 2,
    label = "Numerical: N = {}, h = {:.2e}".format(line_num, 1/(line_num-1)))
plt.plot(x, analytical, linewidth = 2, label = "Analytical: $u(x)=1-(1-e^{-10})x-e^{-10x}$")
plt.xlabel("$x$-coordinate [dim-less]", fontname = "serif", fontsize = 15)
plt.ylabel("$u(x)$ [dim-less]", fontname = "serif", fontsize = 15)
plt.title("Solution to 1D Poisson equation, source: $g(x)=100e^{-10x}$",
    fontname = "serif", fontsize = 16)
plt.legend()
plt.show()

