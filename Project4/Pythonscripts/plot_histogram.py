#PLOT HISTOGRAM 4d)

import numpy as np
import matplotlib.pyplot as mpl
mpl.rcParams.update({'font.size': 20})

#T=2.4
Hist_t2p4 = np.loadtxt('20_4c_hist_t2p4', skiprows=1, dtype=np.float64)
#SigmaE^2 = 1.4063226 * (2.4**2) = 8.10041472


mpl.figure(1)
weights = np.ones_like(Hist_t2p4)/float(len(Hist_t2p4))
mpl.hist(Hist_t2p4, weights=weights, bins=118)
mpl.ylabel('Probability', fontsize=22);
mpl.title('Histogram for Energy at T = 2.4', fontsize=22)
mpl.xlabel('Energy', fontsize=22)
mpl.show()



#T=1.0
Hist_t1_R = np.loadtxt('20_hist_t1R_perN', skiprows=1, dtype=np.float64)
#SigmaE^2 = 0.024205780

mpl.figure(2)
weights = np.ones_like(Hist_t1_R)/float(len(Hist_t1_R))
mpl.hist(Hist_t1_R, weights=weights, bins=6)
mpl.ylabel('Probability' , fontsize=22);
mpl.title('Histogram for Energy at T = 1.0, Random Start Configuration', fontsize=22)
mpl.xlabel('Energy', fontsize=22)
mpl.show()


#T=1.0
Hist_t1_G = np.loadtxt('20_hist_t1G', skiprows=1, dtype=np.float64)
#sigmaE^2 = 0.021664241


mpl.figure(3)
weights = np.ones_like(Hist_t1_G)/float(len(Hist_t1_G))
mpl.hist(Hist_t1_G, weights=weights, bins=6)
mpl.ylabel('Probability', fontsize=22)
mpl.xlabel('Energy', fontsize=22)
mpl.title('Histogram for Energy at T = 1.0, Ground Start Configuration', fontsize=22)
mpl.show()



"""
energy_values=np.array([-800, -792, -788, -784, -780, -776, -772, -768, -764, -760, -756, -752])
number_events=np.array([8.69898e+06, 1.1676e+06, 42488, 80552, 5912, 3879, 407, 152, 12, 6, 2, 1])
total_events=sum(number_events)

#Probability of E=-800 
Prob_min_800 = 8.69898e+06/total_events
print Prob_min_800
"""




