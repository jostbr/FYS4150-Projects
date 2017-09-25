import numpy as np
import matplotlib.pyplot as mpl

n=200

#eig1=[]
#with list I should append

eig_1=np.zeros(n)
eig2=np.zeros(n)
eig3=np.zeros(n)

#print(eig1)

with open("outputfile", "r") as results:
	for line_num, line_string in enumerate(results):
            if (line_num > 0):
                words = line_string.split()
                eig_1[line_num-1] = float(words[0])
                eig2[line_num-1] = float(words[1])
		eig3[line_num-1] = float(words[2])

rho = np.linspace(0, 4.0, n) 

mpl.plot(rho, eig_1, "r-*")
mpl.hold('on')
mpl.plot(rho, eig2, "b-*")
mpl.plot(rho, eig3, "y-*")
mpl.legend(['$\Psi_1$', '$\Psi_2$', '$\Psi_3$'])
mpl.title('Graphical representation of the three first wavevectors')
mpl.xlabel('rho')
mpl.ylabel('Probability; |$\Psi$|^2')
mpl.show()
