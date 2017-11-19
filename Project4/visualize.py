
import numpy as np
import matplotlib.pyplot as plt

# Visualizing results for first part of exercise 4c) (20x20 lattice energy, moment)
# ==================================================================================
filename_01 = "../build-Project4/20x20_T1_RANDOM.txt"
filename_02 = "../build-Project4/20x20_T1_ORDERED.txt"
RANDOM_20x20_T1 = np.loadtxt(filename_01, dtype = np.float64, skiprows = 1)
ORDERED_20x20_T1 = np.loadtxt(filename_02, dtype = np.float64, skiprows = 1)

filename_03 = "../build-Project4/20x20_T24_RANDOM.txt"
filename_04 = "../build-Project4/20x20_T24_ORDERED.txt"
RANDOM_20x20_T24 = np.loadtxt(filename_03, dtype = np.float64, skiprows = 1)
ORDERED_20x20_T24 = np.loadtxt(filename_04, dtype = np.float64, skiprows = 1)

plt.style.use("seaborn")
fig_1, [ax_11, ax_12, ax_13, ax_14] = plt.subplots(figsize = (8, 10), nrows = 4, ncols = 1)

ax_11.plot(RANDOM_20x20_T1[:, 2], RANDOM_20x20_T1[:, 3], label = "Random")
ax_11.plot(ORDERED_20x20_T1[:, 2], ORDERED_20x20_T1[:, 3], label = "Ordered")
ax_11.set_title("Mean energy <$E$> for $T=1.0$", fontsize = 16, fontname = "serif")
#ax_1.set_xlabel("Number of Monte Carlo cycles", fontsize = 14, fontname = "serif")
ax_11.set_ylabel("<$E$>", fontsize = 14, fontname = "serif")
ax_11.legend()

ax_12.plot(RANDOM_20x20_T24[:, 2], RANDOM_20x20_T24[:, 3], label = "Random")
ax_12.plot(ORDERED_20x20_T24[:, 2], ORDERED_20x20_T24[:, 3], label = "Ordered")
ax_12.set_title("Mean energy <$E$> for $T=2.4$", fontsize = 16, fontname = "serif")
#ax_2.set_xlabel("Number of Monte Carlo cycles", fontsize = 14, fontname = "serif")
ax_12.set_ylabel("<$E$>", fontsize = 14, fontname = "serif")
ax_12.legend()

ax_13.plot(RANDOM_20x20_T1[:, 2], RANDOM_20x20_T1[:, 4], label = "Random")
ax_13.plot(ORDERED_20x20_T1[:, 2], ORDERED_20x20_T1[:, 4], label = "Ordered")
ax_13.set_title("Mean magnetization <$|M|$> for $T=1.0$", fontsize = 16, fontname = "serif")
#ax_1.set_xlabel("Number of Monte Carlo cycles", fontsize = 14, fontname = "serif")
ax_13.set_ylabel("<$|M|$>", fontsize = 14, fontname = "serif")
ax_13.legend()

ax_14.plot(RANDOM_20x20_T24[:, 2], RANDOM_20x20_T24[:, 4], label = "Random")
ax_14.plot(ORDERED_20x20_T24[:, 2], ORDERED_20x20_T24[:, 4], label = "Ordered")
ax_14.set_title("Mean magnetization <$|M|$> for $T=2.4$", fontsize = 16, fontname = "serif")
ax_14.set_xlabel("Number of Monte Carlo cycles", fontsize = 14, fontname = "serif")
ax_14.set_ylabel("<$|M|$>", fontsize = 14, fontname = "serif")
ax_14.legend()

fig_1.tight_layout()
# ==================================================================================


# Visualizing results for second part of exercise 4c) (number of accepted moves
# as function of temperature for 20x20 lattice)
# ==================================================================================
filename_04 = "../build-Project4/20x20_TR_RANDOM.txt"
RANDOM_20x20_TR = np.loadtxt(filename_04, dtype = np.float64, skiprows = 1)

plt.style.use("ggplot")
fig_2, ax_21 = plt.subplots(figsize = (7, 7))

ax_21.plot(RANDOM_20x20_TR[:, 1], RANDOM_20x20_TR[:, -1])
ax_21.set_title(r"Metropolis tests passed, $20\times20$ lattice", fontsize = 16, fontname = "serif")
ax_21.set_xlabel("Temperature [J]", fontsize = 14, fontname = "serif")
ax_21.set_ylabel("Total accepted configs", fontsize = 14, fontname = "serif")
fig_2.tight_layout()
# ==================================================================================

# Visualizing rresults for exercise 4d) (probability distribution P(E) for 20x20)
# ==================================================================================
filename_05 = "../build-Project4/energy_dist.txt"
ENERGY_DIST = np.loadtxt(filename_05, dtype = np.float64, skiprows = 1)
#bins, bin_edges = np.histogram(ENERGY_DIST, bins = 100)
#print(hist.shape, len(bin_edges))
print(ENERGY_DIST.shape)

plt.style.use("ggplot")
fig_3, ax_31 = plt.subplots(figsize = (7, 7))

hist, bins, patches = ax_31.hist(ENERGY_DIST, bins = 50, alpha = 0.75, edgecolor = "k")
ax_31.set_title(r"Probability distribution of energy $20\times20$ lattice", fontsize = 16, fontname = "serif")
ax_31.set_xlabel("Energy", fontsize = 14, fontname = "serif")
ax_31.set_ylabel("Number of hits", fontsize = 14, fontname = "serif")
fig_3.tight_layout()
# ==================================================================================





plt.show()