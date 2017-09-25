import numpy as np
import matplotlib.pyplot as mpl

n=100


#eig1=[]
#with list I should append

eig_1=np.zeros(n)
eig2=np.zeros(n)
eig3=np.zeros(n)

#print(eig1)

with open("outputfile.txt", "r") as results:
	for line_num, line_string in enumerate(results):
            if (line_num > 0):
                words = line_string.split()
                eig_1[line_num-1] = float(words[0])
                eig2[line_num-1] = float(words[1])
		eig3[line_num-1] = float(words[2])

rho = np.linspace(0, 10, n) 

mpl.plot(rho, eig_1)
mpl.show()
