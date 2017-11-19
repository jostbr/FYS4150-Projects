import numpy as np
import matplotlib.pyplot as mpl
mpl.rcParams.update({'font.size': 20})

#Random state, T=1.0
MCC = [100, 100200, 300300, 600400, 1000500, 1500600, 2100700, 2800800, 3600900, 4500100, 5501100]
Accepted_Random = [1228, 29873, 87744, 174434, 289672, 435226, 609708, 812095, 1043782, 1304905, 1594354]

#Ground state, T=1.0
#MCC = [100, 100200, 300300, 600400, 1000500, 1500600, 2100700, 2800800, 3600900, 4500100, 5501100]
Accepted_Ground = [46, 29238, 87512, 175151, 290708, 435786, 608897, 811942, 1045514, 1306280, 1596137]

#T=2.4
#MCC = [100, 100200, 300300, 600400, 1000500, 1500600, 2100700, 2800800, 3600900, 4500100, 5501100]
#Ground - Wrong
#Accepted_highT = [10236, 10793467, 32330495, 64777476, 107957303, 167200573, 226614183, 301928455, 388331289, 485379179, 593536817]
#Random
Accepted_highT = [12834, 10877943, 32347804, 64686283, 107635148, 161652528, 226145702, 301271035, 387781863, 484779784, 592658726]

mpl.figure(1)
mpl.subplot(1,2,1)
mpl.plot(MCC, Accepted_Random)
mpl.hold('on')
mpl.xlabel('Total number of Monte Carlo Cycles', fontsize=22)
mpl.ylabel('Accepted New Spin Configurations', fontsize=22)
mpl.plot(MCC, Accepted_Ground)
mpl.legend(['Ground Spin Start Configuration', 'Random Spin Start Configuration'], loc='upper left', fontsize=22)
mpl.title('Accepted configurations, T=1.0, L=20', fontsize=22)
#mpl.show()


mpl.subplot(1,2,2)
mpl.plot(MCC, Accepted_highT)
mpl.hold('on')
mpl.xlabel('Total number of Monte Carlo Cycles', fontsize=22)
mpl.ylabel('Accepted New Spin Configurations', fontsize=22)
#mpl.legend(['Random'], loc='upper center')
mpl.title('Accepted configurations, T=2.4, L=20', fontsize=22)
mpl.show()


Tsweep=np.linspace(1.0, 2.4, 15)
Accepted_random_Tsweep = [291536, 620984, 1187673, 2088411, 3424786, 5356206, 8028852, 11669942, 16537634, 22995396, 31630606, 43289631, 59916869, 83039758, 107658796]
#print Tsweep
Accepted_ground_Tsweep = [289116, 621895, 1188659, 2082875, 3432776, 5350580, 8033567, 11682083, 16539001, 23029236, 31630280, 43417207, 60011156, 83046234, 107608134]

mpl.figure(2)
mpl.plot(Tsweep, Accepted_random_Tsweep)
mpl.hold('on')
mpl.plot(Tsweep, Accepted_ground_Tsweep)
mpl.xlabel('Temperature', fontsize=18)
mpl.ylabel('Accepted New Spin Configurations', fontsize=22)
mpl.legend(['Random Spin Start Configuration', 'Ground Spin Start Configuration'], loc='upper left')
mpl.title('Accepted configurations Vs. Temperature, L=20', fontsize=22)
mpl.axis((1.0, 2.4, 28000, 110000000))
mpl.show()

