
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def animate_wave(t, x, psi):
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
        if (t[i] >= plot_times[0]):
            ax.plot(x, psi[i, :], linewidth = 2, label = "t = {:.2f}".format(t[i]))
            plot_times.pop(0)

    ax.legend()

def compute_analytical_basin(x, t):
    L = 1.0
    k = 4.0*np.pi/L
    omega = -1.0
    psi = np.zeros((t.shape[0], x.shape[0]))

    for i in range(len(t)):
        psi[i, :] = np.sin(k*x)*np.cos(k*x - omega*t[i])

    return psi




filename = "../build-Project5/results.txt"
data = np.loadtxt(filename, dtype = np.float64)

t = data[:, 0]
psi = data[:, 1:]
x = np.linspace(0, 1, psi.shape[1])

psi_analytical = compute_analytical_basin(x, t)
print(psi_analytical.shape)

anim_basin_numerical = animate_wave(t, x, psi)
#anim_basin_analytical = animate_wave(t, x, psi_analytical)
#plot_psi_at_times(x, t, psi, [0, 50, 100, 150])

plt.show()