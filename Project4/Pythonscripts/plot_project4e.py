#Plot Project 4e)
import numpy as np
import matplotlib.pyplot as mpl
mpl.rcParams.update({'font.size': 18})

Spin_40 = np.loadtxt('40_4e_new', skiprows=1, dtype=np.float64)
#Temperature    Total MCC       E            Cv             M            Chi            |M| 

mpl.figure(1)
mpl.subplot(4,1,1)
mpl.plot(Spin_40[:, 0], Spin_40[:, 2], '-b')
mpl.hold('on')
mpl.ylabel('Energy', fontsize=18)
mpl.title('Expectation Values per Spin for a 40x40 spin lattice', fontsize=20)

mpl.subplot(4,1,2)
mpl.plot(Spin_40[:, 0], Spin_40[:, 6], '-b')
mpl.ylabel('Mean Magnetization', fontsize=18)

mpl.subplot(4,1,3)
mpl.plot(Spin_40[:, 0], Spin_40[:, 3], '-b')
mpl.ylabel('Heat Capacity', fontsize=18)

mpl.subplot(4,1,4)
mpl.plot(Spin_40[:, 0], Spin_40[:, 5], '-b')
mpl.xlabel('Temperature', fontsize=18)
mpl.ylabel('Susceptibility', fontsize=18)

mpl.show()



Spin_60 = np.loadtxt('60_4e_new_new', skiprows=1, dtype=np.float64)
#Temperature    Total MCC       E            Cv             M            Chi            |M| 

mpl.figure(2)
mpl.subplot(4,1,1)
mpl.plot(Spin_60[:, 0], Spin_60[:, 2], '-b')
mpl.hold('on')
mpl.ylabel('Energy', fontsize=18)
mpl.title('Expectation Values per Spin for a 60x60 spin lattice', fontsize=20)

mpl.subplot(4,1,2)
mpl.plot(Spin_60[:, 0], Spin_60[:, 6], '-b')
mpl.ylabel('Mean Magnetization', fontsize=18)

mpl.subplot(4,1,3)
mpl.plot(Spin_60[:, 0], Spin_60[:, 3], '-b')
mpl.ylabel('Heat Capacity', fontsize=18)

mpl.subplot(4,1,4)
mpl.plot(Spin_60[:, 0], Spin_60[:, 5], '-b')
mpl.xlabel('Temperature', fontsize=18)
mpl.ylabel('Susceptibility', fontsize=18)

mpl.show()



Spin_80 = np.loadtxt('80_4e', skiprows=1, dtype=np.float64)
#Temperature    Total MCC       E            Cv             M            Chi            |M| 

mpl.figure(1)
mpl.subplot(4,1,1)
mpl.plot(Spin_80[:, 0], Spin_80[:, 2], '-b')
mpl.hold('on')
mpl.ylabel('Energy', fontsize=18)
mpl.title('Expectation Values per Spin for a 80x80 spin lattice', fontsize=20)

mpl.subplot(4,1,2)
mpl.plot(Spin_80[:, 0], Spin_80[:, 6], '-b')
mpl.ylabel('Mean Magnetization', fontsize=18)

mpl.subplot(4,1,3)
mpl.plot(Spin_80[:, 0], Spin_80[:, 3], '-b')
mpl.ylabel('Heat Capacity', fontsize=18)

mpl.subplot(4,1,4)
mpl.plot(Spin_80[:, 0], Spin_80[:, 5], '-b')
mpl.xlabel('Temperature', fontsize=18)
mpl.ylabel('Susceptibility', fontsize=18)

mpl.show()




Spin_100 = np.loadtxt('100_4e_new', skiprows=1, dtype=np.float64)
#Temperature    Total MCC       E            Cv             M            Chi            |M| 

mpl.figure(1)
mpl.subplot(4,1,1)
mpl.plot(Spin_100[:, 0], Spin_100[:, 2], '-b')
mpl.hold('on')
mpl.ylabel('Energy', fontsize=18)
mpl.title('Expectation Values per Spin for a 100x100 spin lattice', fontsize=20)

mpl.subplot(4,1,2)
mpl.plot(Spin_100[:, 0], Spin_100[:, 6], '-b')
mpl.ylabel('Mean Magnetization', fontsize=18)

mpl.subplot(4,1,3)
mpl.plot(Spin_100[:, 0], Spin_100[:, 3], '-b')
mpl.ylabel('Heat Capacity', fontsize=18)

mpl.subplot(4,1,4)
mpl.plot(Spin_100[:, 0], Spin_100[:, 5], '-b')
mpl.xlabel('Temperature', fontsize=18)
mpl.ylabel('Susceptibility', fontsize=18)

mpl.show()

