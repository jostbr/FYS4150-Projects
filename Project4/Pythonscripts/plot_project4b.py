import numpy as np
import matplotlib.pyplot as mpl
mpl.rcParams.update({'font.size': 20})

#Header of .txt file is 
#Temperature   Tot_MCC    E      Cv       M     Chi      |M| 

Spin_2R = np.loadtxt('2_T1R', skiprows=1, dtype=np.float64)
Spin_2G = np.loadtxt('2_T1G', skiprows=1, dtype=np.float64)


mpl.figure(1)
Analytical_E_value = -1.99598
Analytical_E = np.ones((len(Spin_2G[:, 1]),1))
Analytical_E = Analytical_E*Analytical_E_value


mpl.figure(1)
mpl.plot(Spin_2G[:, 1], Spin_2G[:, 2], '-b')
mpl.hold('on')
mpl.plot(Spin_2R[:, 1], Spin_2R[:, 2], '-g')
mpl.plot(Spin_2G[:, 1], Analytical_E, '-r')
mpl.xlabel('Number of Monte Carlo Cycles', fontsize=22)
mpl.ylabel('Expectation Value for Energy per Spin', fontsize=22)
mpl.legend(['Numerical Result Ground Start Configuration of SpinLattice','Numerical Result Random Start Configuration of SpinLattice','Analytical Result'], fontsize=20)
mpl.axis([0, 10000000, -1.997, -1.9945])
mpl.title('2x2 spin lattice, T=1.0', fontsize=22)
mpl.show()

##################################################################################

mpl.figure(2)
Analytical_M_value =  0.998661
Analytical_M = np.ones(len(Spin_2G[:, 1]))
Analytical_M = Analytical_M_value*Analytical_M

mpl.plot(Spin_2G[:, 1], Spin_2G[:, 6], '-b')
mpl.hold('on')
mpl.plot(Spin_2R[:, 1], Spin_2R[:, 6], '-g')
mpl.plot(Spin_2G[:, 1], Analytical_M, '-r')
mpl.xlabel('Number of Monte Carlo Cycles', fontsize=22)
mpl.ylabel('Expectation Value for Absolut Magnetization per Spin', fontsize=22)
mpl.legend(['Numerical Result Ground Start Configuration of SpinLattice','Numerical Result Random Start Configuration of SpinLattice','Analytical Result'], fontsize=20)
mpl.axis([0, 10000000, 0.9982, 0.9990])
mpl.title('2x2 spin lattice, T=1.0', fontsize=22)
mpl.show()


##########################################################################


Analytical_Cv_value =  0.0320823
Analytical_Cv = np.ones(len(Spin_2G[:, 1]))
Analytical_Cv = Analytical_Cv_value*Analytical_Cv

mpl.figure(3)
mpl.plot(Spin_2G[:, 1], Spin_2G[:, 3], '-b')
mpl.hold('on')
mpl.plot(Spin_2R[:, 1], Spin_2R[:, 3], '-g')
mpl.plot(Spin_2G[:, 1], Analytical_Cv, '-r')
mpl.xlabel('Number of Monte Carlo Cycles', fontsize=22)
mpl.ylabel('Expectation Value for Heat Capacity per Spin', fontsize=22)
mpl.legend(['Numerical Result Ground Start Configuration of SpinLattice','Numerical Result Random Start Configuration of SpinLattice','Analytical Result'], fontsize=20)
mpl.axis([0, 10000000, 0.025, 0.045])
mpl.title('2x2 spin lattice, T=1.0', fontsize=22)
mpl.show()

###################################################################################

Analytical_Chi_value = 0.00401074
Analytical_Chi = np.ones(len(Spin_2G[:, 1]))
Analytical_Chi = Analytical_Chi_value*Analytical_Chi


mpl.figure(4)
mpl.plot(Spin_2G[:, 1], Spin_2G[:, 5], '-b')
mpl.hold('on')
mpl.plot(Spin_2R[:, 1], Spin_2R[:, 5], '-g')
mpl.plot(Spin_2G[:, 1], Analytical_Chi, '-r')
mpl.xlabel('Number of Monte Carlo Cycles', fontsize=22)
mpl.ylabel('Expectation Value for susceptibility per Spin', fontsize=22)
mpl.legend(['Numerical Result Ground Start Configuration of SpinLattice','Numerical Result Random Start Configuration of SpinLattice','Analytical Result'], fontsize=20)
mpl.axis([0, 10000000, 0.003, 0.005])
mpl.title('2x2 spin lattice, T=1.0', fontsize=20)
mpl.show()

