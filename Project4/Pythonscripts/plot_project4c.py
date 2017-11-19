import numpy as np
import matplotlib.pyplot as mpl
mpl.rcParams.update({'font.size': 20})

#Header of .txt file is 
#Temperature   Tot_MCC    E      Cv       M     Chi      |M| 

#Load the datasets
Spin_20R_T1 = np.loadtxt('20_T1_Random_4c', skiprows=1, dtype=np.float64)	    
Spin_20G_T1 = np.loadtxt('20_4c_t1G_final', skiprows=1, dtype=np.float64) 

#Spin_20_T2p4 = np.loadtxt('20_4c_t2p4_new', skiprows=1, dtype=np.float64)	  
Spin_20_T2p4 = np.loadtxt('20_T2p4_4c', skiprows=1, dtype=np.float64)

#####################################################################

Spin_20G2_T1 = np.loadtxt('20_T1G22_4c', skiprows=1, dtype=np.float64) 
Spin_20R2_T1 = np.loadtxt('20_T1R2_4c', skiprows=1, dtype=np.float64)

mpl.figure(4)
mpl.plot(Spin_20R2_T1[:,1], Spin_20R2_T1[:,2], '-b')
mpl.hold('on')
mpl.plot(Spin_20G2_T1[:,1], Spin_20G2_T1[:,2], '-g')
mpl.xlabel('Number of Monte Carlo Cycles', fontsize=20)
mpl.ylabel('Mean Energy', fontsize=20)
mpl.legend(['Random Start Configuration','Ground Start Configuration'])
mpl.title('Equlibriation time, L=20, T=1.0', fontsize=20)
mpl.axis([0, 16000, -2.0, -1.99])
mpl.show()


"""

mpl.figure(3)
mpl.subplot(4,1,1)
mpl.plot(Spin_20_T2p4[:, 1], Spin_20_T2p4[:, 2], '-b')
mpl.ylabel('Mean Energy', fontsize=18)
mpl.title('Expectation Values per Spin for a 20x20 spin lattice, T=2.4', fontsize=20)
mpl.axis([0, 20400000, -1.24, -1.2])
mpl.hold('on')
mpl.grid('on')

mpl.subplot(4,1,2)
mpl.plot(Spin_20_T2p4[:, 1], Spin_20_T2p4[:, 6], '-b')
mpl.ylabel('Mean Magnetization', fontsize=18)
mpl.axis([0, 20400000, 0.40, 0.48])
mpl.grid('on')

mpl.subplot(4,1,3)
mpl.plot(Spin_20_T2p4[:, 1], Spin_20_T2p4[:, 3], '-b')
mpl.ylabel('Heat Capacity', fontsize=18)
mpl.axis([0,20400000, 1.3, 1.5])
mpl.grid('on')

mpl.subplot(4,1,4)
mpl.plot(Spin_20_T2p4[:, 1], Spin_20_T2p4[:, 5], '-b')
mpl.xlabel('Number of Monte Carlo Cycles', fontsize=18)
mpl.ylabel('Susceptibility', fontsize=18)
mpl.axis([0, 20400000, 8.6, 9.2])
mpl.grid('on')
mpl.show()

#######################################################################

mpl.figure(1)
mpl.subplot(4,1,1)
mpl.plot(Spin_20G_T1[:, 1], Spin_20G_T1[:, 2], '-b')
mpl.hold('on')
mpl.grid('on')
mpl.axis([0, 10000000, -1.998, -1.997])
mpl.ylabel('Mean Energy', fontsize=18)
mpl.title('Expectation Values per Spin for a 20x20 spin lattice, Ground Start Configuration, T=1.0', fontsize=20)

mpl.subplot(4,1,2)
mpl.plot(Spin_20G_T1[:, 1], Spin_20G_T1[:, 6], '-b')
mpl.grid('on')
mpl.axis([0, 10000000, 0.9992, 0.99935 ])

mpl.ylabel('Mean Magnetization', fontsize=18)


mpl.subplot(4,1,3)
mpl.plot(Spin_20G_T1[:, 1], Spin_20G_T1[:, 3], '-b')
mpl.grid('on')
mpl.ylabel('Heat Capacity', fontsize=18)
mpl.axis([0, 10000000, 0.015, 0.0275])


mpl.subplot(4,1,4)
mpl.plot(Spin_20G_T1[:, 1], Spin_20G_T1[:, 5], '-b')
mpl.grid('on')
mpl.xlabel('Number of Monte Carlo Cycles', fontsize=18)
mpl.ylabel('Susceptibility', fontsize=18)
mpl.axis([0, 10000000, 0.0, 0.002])
mpl.show()


#######################################################################


mpl.figure(2)

mpl.subplot(4,1,1)
mpl.hold('on')
mpl.grid('on')
mpl.plot(Spin_20R_T1[:, 1], Spin_20R_T1[:, 2], '-g')
mpl.ylabel('Mean Energy', fontsize=18)
mpl.title('Expectation Values per Spin for a 20x20 spin lattice, Random Start Configuration, T=1.0', fontsize=20)
mpl.axis([0, 20000000, -1.998, -1.997])

mpl.subplot(4,1,2)
mpl.plot(Spin_20R_T1[:, 1], Spin_20R_T1[:, 6], '-g')
mpl.grid('on')
mpl.ylabel('Mean Magnetization', fontsize=18)
mpl.axis([0, 20000000, 0.99925, 0.99935 ])

mpl.subplot(4,1,3)
mpl.grid('on')
mpl.plot(Spin_20R_T1[:, 1], Spin_20R_T1[:, 3], '-g')
mpl.ylabel('Heat Capacity', fontsize=18)
mpl.axis([0, 20000000, 0.018, 0.025])

mpl.subplot(4,1,4)
mpl.grid('on')
mpl.plot(Spin_20R_T1[:, 1], Spin_20R_T1[:, 5], '-g')
mpl.xlabel('Number of Monte Carlo Cycles', fontsize=18)
mpl.ylabel('Susceptibility', fontsize=18)
mpl.axis([0, 20000000, -0.002, 0.002])
mpl.show()

"""






