
import numpy as np
import matplotlib.pyplot as plt

#with open("../build-Project3/results.txt", "r") as file_result:
#	for line_num, linestring in enumerate(file_result):

filename = "../build-Project3/results.txt"
data = np.loadtxt(filename, dtype = np.float64, skiprows = 1)
print(data.shape)
plt.plot(data[:, 1], data[:, 2])
plt.show()
