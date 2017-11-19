
import numpy as np
import matplotlib.pyplot as plt

# Visualizing results for first part of exercise 4c) (20x20 lattice energy, moment)
# ==================================================================================
filename_01 = "benchmarks/20x20_T1_RANDOM.txt"
filename_02 = "benchmarks/20x20_T1_ORDERED.txt"
RANDOM_20x20_T1 = np.loadtxt(filename_01, dtype = np.float64, skiprows = 1)
ORDERED_20x20_T1 = np.loadtxt(filename_02, dtype = np.float64, skiprows = 1)

filename_03 = "benchmarks/20x20_T24_RANDOM.txt"
filename_04 = "benchmarks/20x20_T24_ORDERED.txt"
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
filename_04 = "benchmarks/20x20_TR_RANDOM.txt"
RANDOM_20x20_TR = np.loadtxt(filename_04, dtype = np.float64, skiprows = 1)

plt.style.use("ggplot")
fig_2, ax_21 = plt.subplots(figsize = (7, 7))

ax_21.plot(RANDOM_20x20_TR[:, 1], RANDOM_20x20_TR[:, -1])
ax_21.set_title(r"Metropolis tests passed, $20\times20$ lattice", fontsize = 16, fontname = "serif")
ax_21.set_xlabel("Temperature [J]", fontsize = 14, fontname = "serif")
ax_21.set_ylabel("Total accepted configs", fontsize = 14, fontname = "serif")

fig_2.tight_layout()
# ==================================================================================

# Visualizing results for exercise 4d) (probability distribution P(E) for 20x20)
# ==================================================================================
filename_05 = "benchmarks/20x20_T1_RANDOM_ENERGY_DIST.txt"
filename_06 = "benchmarks/20x20_T24_RANDOM_ENERGY_DIST.txt"
ENERGY_DIST_T1 = np.loadtxt(filename_05, dtype = np.float64, skiprows = 0)
ENERGY_DIST_T24 = np.loadtxt(filename_06, dtype = np.float64, skiprows = 0)
hist_T1, bin_edges_T1 = np.histogram(ENERGY_DIST_T1, bins = 7, density = True)
hist_T24, bin_edges_T24 = np.histogram(ENERGY_DIST_T24, bins = 18, density = True)

plt.style.use("seaborn-dark")
fig_3, [ax_31, ax_32] = plt.subplots(figsize = (11, 7), nrows = 1, ncols = 2)

#hist, bins, patches = ax_31.hist(ENERGY_DIST, bins = 6, align = "left", alpha = 0.75, edgecolor = "k")
ax_31.bar(bin_edges_T1[:-1], hist_T1, width = np.diff(bin_edges_T1)[0], align = "center", edgecolor = "k", alpha = 0.75)
ax_31.set_title(r"Probability distribution $P(E)$ for $T=1.0$", fontsize = 15, fontname = "serif")
ax_31.set_xlabel("Energy per spin", fontsize = 12, fontname = "serif")
ax_31.set_ylabel("Normilized number of occurences", fontsize = 12, fontname = "serif")

ax_32.bar(bin_edges_T24[:-1], hist_T24, width = np.diff(bin_edges_T24)[0], align = "center", edgecolor = "k", alpha = 0.75)
ax_32.set_title(r"Probability distribution $P(E)$ for $T=2.4$", fontsize = 15, fontname = "serif")
ax_32.set_xlabel("Energy per spin", fontsize = 12, fontname = "serif")
ax_32.set_ylabel("Normilized number of occurences", fontsize = 12, fontname = "serif")
ax_31.grid(True)
ax_32.grid(True)

fig_3.tight_layout()
# ==================================================================================

# Visualizing results for exercise 4e) (phase transitions)
# ==================================================================================
filename_07 = "benchmarks/PR_20x20_TR.txt"
filename_08 = "benchmarks/PR_40x40_TR.txt"
filename_09 = "benchmarks/PR_60x60_TR.txt"
filename_10 = "benchmarks/PR_80x80_TR.txt"
PR_20x20 = np.loadtxt(filename_07, dtype = np.float64, skiprows = 1)
PR_40x40 = np.loadtxt(filename_08, dtype = np.float64, skiprows = 1)
PR_60x60 = np.loadtxt(filename_09, dtype = np.float64, skiprows = 1)
PR_80x80 = np.loadtxt(filename_10, dtype = np.float64, skiprows = 1)

plt.style.use("ggplot")
fig_4, [ax_41, ax_42, ax_43, ax_44] = plt.subplots(figsize = (8, 10), nrows = 4, ncols = 1)

ax_41.plot(PR_20x20[:, 1], PR_20x20[:, 3], label = r"$20\times20$")
ax_41.plot(PR_40x40[:, 1], PR_40x40[:, 3], label = r"$40\times40$")
ax_41.plot(PR_60x60[:, 1], PR_60x60[:, 3], label = r"$60\times60$")
ax_41.plot(PR_80x80[:, 1], PR_80x80[:, 3], label = r"$80\times80$")
ax_41.set_title("Mean energy for various lattices", fontsize = 16, fontname = "serif")
#ax_41.set_xlabel("Temperature [J]", fontsize = 14, fontname = "serif")
ax_41.set_ylabel("<$E$>", fontsize = 14, fontname = "serif")
ax_41.legend()

ax_42.plot(PR_20x20[:, 1], PR_20x20[:, 4], label = r"$20\times20$")
ax_42.plot(PR_40x40[:, 1], PR_40x40[:, 4], label = r"$40\times40$")
ax_42.plot(PR_60x60[:, 1], PR_60x60[:, 4], label = r"$60\times60$")
ax_42.plot(PR_80x80[:, 1], PR_80x80[:, 4], label = r"$80\times80$")
ax_42.set_title("Mean magnetization for various lattices", fontsize = 16, fontname = "serif")
#ax_42.set_xlabel("Temperature [J]", fontsize = 14, fontname = "serif")
ax_42.set_ylabel("<$|M|$>", fontsize = 14, fontname = "serif")
ax_42.legend()

ax_43.plot(PR_20x20[:, 1], PR_20x20[:, 8], label = r"$20\times20$")
ax_43.plot(PR_40x40[:, 1], PR_40x40[:, 8], label = r"$40\times40$")
ax_43.plot(PR_60x60[:, 1], PR_60x60[:, 8], label = r"$60\times60$")
ax_43.plot(PR_80x80[:, 1], PR_80x80[:, 8], label = r"$80\times80$")
ax_43.set_title("Specific heat for various lattices", fontsize = 16, fontname = "serif")
#ax_43.set_xlabel("Temperature [J]", fontsize = 14, fontname = "serif")
ax_43.set_ylabel("$C_v$", fontsize = 14, fontname = "serif")
ax_43.legend()

ax_44.plot(PR_20x20[:, 1], PR_20x20[:, 9], label = r"$20\times20$")
ax_44.plot(PR_40x40[:, 1], PR_40x40[:, 9], label = r"$40\times40$")
ax_44.plot(PR_60x60[:, 1], PR_60x60[:, 9], label = r"$60\times60$")
ax_44.plot(PR_80x80[:, 1], PR_80x80[:, 9], label = r"$80\times80$")
ax_44.set_title("Magentic susceptibility for various lattices", fontsize = 16, fontname = "serif")
ax_44.set_xlabel("Temperature [J]", fontsize = 14, fontname = "serif")
ax_44.set_ylabel("$\chi$", fontsize = 14, fontname = "serif")
ax_44.legend()

fig_4.tight_layout()
# ==================================================================================

plt.show()    # Show everything