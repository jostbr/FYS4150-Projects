import numpy as np
import matplotlib.pyplot as mpl
#Stuff that need to be manually changed for each plot
n=500
omega=1.0

m_e=9.10938356e-31	#kg
hbar=1.0545718e-34	#Js


k=m_e*omega*omega
alpha4= ((hbar*hbar)/m_e*k)
alpha = alpha4**0.25

rho = np.linspace(0, 60.0, n)
r=rho*alpha 



eig1=np.zeros(n)
eig2=np.zeros(n)
eig3=np.zeros(n)

#print(eig1)

with open("outputfile", "r") as results:
	for line_num, line_string in enumerate(results):
	    if (line_num > 0):
                words = line_string.split()
                eig1[line_num-1] = float(words[0])
                eig2[line_num-1] = float(words[1])
		eig3[line_num-1] = float(words[2])





mpl.plot(r, eig1, "r-*")
mpl.hold('on')
mpl.plot(r, eig2, "b-*")
mpl.plot(r, eig3, "y-*")
mpl.legend(['$\Psi_1$', '$\Psi_2$', '$\Psi_3$'])
mpl.title('Graphical representation of the three first wavevectors')
mpl.xlabel('rho')
mpl.ylabel('Probability; |$\Psi$|^2')
mpl.show()
