
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def get_data(filename):
    data = np.loadtxt(filename, dtype = np.float64)

    t = data[:, 0]
    psi = data[:, 1:]
    x = np.linspace(0, 1, psi.shape[1])

    return x, t, psi

def animate_wave(x, t, psi):
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (8, 8))
    ax = plt.axes(xlim = (x[0], x[-1]), ylim = (psi.min(), psi.max()))
    ax.set_title("Streamfunction $\psi$ after {:.2f} time".format(t[0]),
        fontname = "serif", fontsize = 18)
    ax.set_xlabel("x [dim-less]", fontname = "serif", fontsize = 13)
    ax.set_ylabel("$\psi(x,t)$ [dim-less]", fontname = "serif", fontsize = 13)

    wave = ax.plot(x, psi[0, :], linewidth = 2)[0]


    def animate(frame):
        wave.set_data(x, psi[frame, :])
        ax.set_title("Streamfunction $\psi(x,t)$ after {:.2f} time".format(t[frame]),
            fontname = "serif", fontsize = 18)
        return []

    anim = animation.FuncAnimation(fig, animate, frames = psi.shape[0], interval = 10, blit = True)
    return anim

def plot_psi_at_times(x, t, psi, plot_times):
    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize = (8, 8))
    ax.set_title("Streamfunction $\psi(x,t)$ at selected times", fontname = "serif", fontsize = 18)
    ax.set_xlabel("x [dim-less]", fontname = "serif", fontsize = 13)
    ax.set_ylabel("$\psi(x,t)$ [dim-less]", fontname = "serif", fontsize = 13)

    for i in range(len(t)):
        if (plot_times and t[i] >= plot_times[0]):
            ax.plot(x, psi[i, :], linewidth = 2, label = "t = {:.2f}".format(t[i]))
            plot_times.pop(0)

    ax.legend()

def compute_analytical_basin(x, t):
    L = 1.0
    k = 2.0*np.pi/L
    omega = -0.5
    psi = np.zeros((t.shape[0], x.shape[0]))

    for i in range(len(t)):
        psi[i, :] = np.sin(k*x)*np.cos(k*x - omega*t[i])

    return psi

def compute_analytical_periodic(x, t):
    L = 1.0
    k = 2.0*np.pi/L
    omega = -1.0
    psi = np.zeros((t.shape[0], x.shape[0]))

    for i in range(len(t)):
        psi[i, :] = np.cos(k*x - omega*t[i])

    return psi


# Plotting and animating basin euler solution
# ====================================================================================
filename_beuler = "../build-Project5/sine_basin_euler.txt"
x_beuler, t_beuler, psi_beuler = get_data(filename_beuler)
#psi_basin_analytical = compute_analytical_basin(x_beuler, t_beuler)

anim_beuler_numerical = animate_wave(x_beuler, t_beuler, psi_beuler)
#anim_basin_analytical = animate_wave(t, x, psi_analytical)
#plot_psi_at_times(x_basin, t_basin, psi_basin, [0, 50, 100, 150])
# ====================================================================================

plt.show()

# Plotting and animating basin leapfrog solution
# ====================================================================================
filename_bleap = "../build-Project5/sine_basin_leapfrog.txt"
x_bleap, t_bleap, psi_bleap = get_data(filename_bleap)
#psi_basin_analytical = compute_analytical_basin(x_bleap, t_bleap)

anim_bleap_numerical = animate_wave(x_bleap, t_bleap, psi_bleap)
#anim_basin_analytical = animate_wave(t, x, psi_analytical)
#plot_psi_at_times(x_basin, t_basin, psi_basin, [0, 50, 100, 150])
# ====================================================================================

plt.show()

# Plotting and animating periodic euler solution
# ====================================================================================
filename_peuler = "../build-Project5/sine_periodic_euler.txt"
x_peuler, t_peuler, psi_peuler = get_data(filename_peuler)

anim_peuler_numerical = animate_wave(x_peuler, t_peuler, psi_peuler)
#plot_psi_at_times(x_periodic, t_periodic, psi_periodic, [0, 5, 10, 50, 100])
# ====================================================================================

plt.show()

# Plotting and animating periodic leapfrog solution
# ====================================================================================
filename_pleap = "../build-Project5/sine_periodic_leapfrog.txt"
x_pleap, t_pleap, psi_pleap = get_data(filename_pleap)

anim_pleap_numerical = animate_wave(x_pleap, t_pleap, psi_pleap)
#plot_psi_at_times(x_periodic, t_periodic, psi_periodic, [0, 5, 10, 50, 100])
# ====================================================================================

plt.show()

psi_basin_analytical = compute_analytical_basin(x_bleap, t_bleap)
anim_basin_analytical = animate_wave(x_bleap, t_bleap, psi_basin_analytical)

psi_periodic_analytical = compute_analytical_periodic(x_pleap, t_pleap)
anim_periodic_analytical = animate_wave(x_pleap, t_pleap, psi_periodic_analytical)

plt.show()
