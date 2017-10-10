
import numpy as np
import matplotlib.pyplot as plt

#with open("../build-Project3/results.txt", "r") as file_result:
#	for line_num, linestring in enumerate(file_result):

filename = "../build-Project3/{}.txt"

data_earth = np.loadtxt(filename.format("earth"), dtype = np.float64, skiprows = 1)
#data_jupiter = np.loadtxt(filename.format("jupiter"), dtype = np.float64, skiprows = 1)
print(data_earth.shape)

plt.plot(data_earth[:, 1], data_earth[:, 2])
#plt.plot(data_jupiter[:, 1], data_jupiter[:, 2])
plt.show()
