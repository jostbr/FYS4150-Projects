
# Script for plotting time usage of the general and
# special algorithm up to n = 10^6 for the grid points.
# Depends on having generated timing results through the
# benchamrk() function in tridiag.cpp.

import numpy as np
import matplotlib.pyplot as plt

def get_timing_data(filename):
    N = list()
    time = list()

    with open(filename, "r") as timing_file:
        for line_num, line_string in enumerate(timing_file):
            if (line_num > 0):
                words = line_string.split()
                N.append(float(words[0]))
                time.append(float(words[1]))

    return np.array(N), np.array(time)

if (__name__ == "__main__"):
    path = "../build-Project1/{}"
    N_gen, time_gen = get_timing_data(path.format("timing_general.txt"))
    N_spe, time_spe = get_timing_data(path.format("timing_special.txt"))

    plt.figure(figsize = (8, 7))
    plt.style.use("seaborn")
    plt.plot(np.array(N_gen)[:-3], np.array(time_gen)[:-3], linewidth = 2, label = "General")
    plt.plot(np.array(N_spe)[:-3], np.array(time_spe)[:-3], linewidth = 2, label = "Special")
    plt.title("Time usage for the two algorithms", fontname = "serif", fontsize = 17)
    plt.xlabel("Size of matrix (NxN)", fontname = "serif", fontsize = 14)
    plt.ylabel("Runtime of algorithm [s]", fontname = "serif", fontsize = 14)
    plt.legend()


    ratio = time_gen/time_spe

    plt.figure(figsize = (8, 7))
    plt.style.use("seaborn")
    plt.plot(N_gen, ratio, linewidth = 2)
    plt.title("Ratio of time usage for the two algorithms", fontname = "serif", fontsize = 18)
    plt.xlabel("Size of matrix (NxN)", fontname = "serif", fontsize = 16)
    plt.ylabel("Ratio of runtimes", fontname = "serif", fontsize = 16)
    plt.show()


