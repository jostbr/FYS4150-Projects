
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def animate_wave(t, x, psi):
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (8, 8))
    ax = plt.axes(xlim = (x[0], x[-1]))
    ax.set_title("Streamfunction $\psi$ after {:.2f} time".format(t[0]),
        fontname = "serif", fontsize = 18)
    ax.set_xlabel("x [dim-less]", fontname = "serif", fontsize = 12)
    ax.set_ylabel("y [dim-less]", fontname = "serif", fontsize = 12)

    wave = ax.plot(x, psi[0, :], linewidth = 2)[0]


    def animate(frame):
        wave.set_data(x, psi[frame, :])
        ax.set_title("Streamfunction $\psi$ after {:.2f} time".format(t[frame]),
            fontname = "serif", fontsize = 18)
        return []

    anim = animation.FuncAnimation(fig, animate, frames = psi.shape[0], interval = 10, blit = True)
    return anim

filename = "../build-Project5/results.txt"
data = np.loadtxt(filename, dtype = np.float64)

t = data[:, 0]
psi = data[:, 1:]

x = np.linspace(0, 1, psi.shape[1])

#plt.plot(x, psi[0, :])

anim = animate_wave(t, x, psi)
plt.show()