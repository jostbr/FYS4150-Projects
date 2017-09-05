
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

upper_exponent = int(sys.argv[1])
problem_sizes = np.logspace(1, upper_exponent, upper_exponent,
    endpoint = True).astype(int)
num_boundary_points = 2

os.system("g++ -o tridiag tridiag.cpp")
os.system("./tridiag {}".format(upper_exponent))

plt.figure(figsize = (8, 7))
plt.style.use("seaborn")

for size in problem_sizes:
    full_system_size = size + num_boundary_points     # Need to account for boundary points as well
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

        plt.plot(x, numerical, linewidth = 2, label = "numerical {}".format(size))

plt.plot(x, analytical, linewidth = 2, label = "analyical")
plt.xlabel("$x$-coordinate", fontname = "serif", fontsize = 16)
plt.ylabel("$u(x)$", fontname = "serif", fontsize = 16)
plt.title("Solution to 1D Poisson eq. (gen. alg.)",
    fontname = "serif", fontsize = 18)
plt.legend()
plt.show()
