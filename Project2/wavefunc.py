
import numpy as np
import matplotlib.pyplot as mpl

# Stuff that need to be manually changed for each plot
n=202
omega=5.0

#m_e=9.10938356e-31	#kg
#hbar=1.0545718e-34	#Js = 
#hbar=6.582119514e-16 	#eV s

#k=m_e*omega*omega

#beta=1.44e-9	#eV m

#alpha= ((hbar*hbar)/m_e*beta)
#alpha = alpha4**0.25

rho = np.linspace(0, 2.0, n)
#r=rho*alpha 

# To hold the three eigenvectors.
eig1=np.zeros(n)
eig2=np.zeros(n)
eig3=np.zeros(n)

#print(eig1)

with open("../build-Project2/results_int_omega5_200.txt", "r") as results:
    for line_num, line_string in enumerate(results):
        if (line_num > 0):
                words = line_string.split()
                eig1[line_num-1] = float(words[0])
                eig2[line_num-1] = float(words[1])
                eig3[line_num-1] = float(words[2])


mpl.figure(figsize = (7, 6))
mpl.plot(rho, eig1, linewidth = 2)
mpl.plot(rho, eig2, linewidth = 2)
mpl.plot(rho, eig3, linewidth = 2)
mpl.legend(['$\Psi_1$', '$\Psi_2$', '$\Psi_3$'])
mpl.title('Plot of the first three wavevectors, $\omega_r={}$'.format(omega),
    fontname = "serif", fontsize = 16)
mpl.xlabel('rho [dim-less]', fontname = "serif", fontsize = 13)
mpl.ylabel('Probability; |$\Psi$|^2', fontname = "serif", fontsize = 13)
mpl.show()
