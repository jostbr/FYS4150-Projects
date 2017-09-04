
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

os.system("g++ -o tridiag tridiag.cpp")
os.system("./tridiag {}".format(sys.argv[1]))

x = list()
numerical = list()
analytical = list()
rel_error = list()

with open("output.txt", "r") as results:
    for line_num, line_string in enumerate(results):
        words = line_string.split()
        x.append(float(words[0]))
        numerical.append(float(words[1]))
        analytical.append(float(words[2]))
        rel_error.append(float(words[3]))

plt.style.use("ggplot")
plt.plot(x, numerical, linewidth = 2, label = "numerical")
plt.plot(x, analytical, linewidth = 2, label = "analyical")
plt.legend()
plt.show()
