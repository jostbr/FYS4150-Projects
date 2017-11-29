
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

def get_data(filename):
    data = np.loadtxt(filename, dtype = np.float64)

    t = data[:, 0]
    psi = data[:, 1:]
    x = np.linspace(0, 1, 41)
    y = np.linspace(0, 1, 41)

    return x, y, t, psi

def contour_results(x, y, psi):
	fig = plt.figure(figsize = (9, 7))
	plt.title("Contour field of $\psi$", fontname = "serif", fontsize = 17)
	plt.xlabel("x [dim-less]", fontname = "serif", fontsize = 12)
	plt.ylabel("y [dim-less]", fontname = "serif", fontsize = 12)
	CS = plt.contourf(x, y, psi.transpose(), 20, cmap = plt.cm.RdBu_r)
	plt.colorbar(CS, orientation = "vertical")

def animate_pcolormesh(x, y, t, psi, filename_save = None):
    """Function that generates an animation of psi over time."""
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (8, 6))
    ax = fig.add_subplot(1, 1, 1)
    
    plot_title = "Streamfunction $\psi(x,y,t)$ at time {}"
    ax.set_title(plot_title.format(t[0]), fontname = "serif", fontsize = 17)
    ax.set_xlabel("x [dim-less]", fontname = "serif", fontsize = 12)
    ax.set_ylabel("y [dim-less]", fontname = "serif", fontsize = 12)
    pmesh = ax.pcolormesh(x, y, psi[0,:-1,:-1], vmin = -1, vmax = 1, cmap = plt.cm.magma)
    fig.colorbar(pmesh, orientation = "vertical")

    def update_frame(frame):
        """Function for updating frame in animation."""
        ax.set_title(plot_title.format(t[frame]),
            fontname = "serif", fontsize = 16)
        pmesh.set_array(psi[frame,:-1,:-1].flatten())     # Set new image
        return pmesh,

    anim = animation.FuncAnimation(fig, update_frame,
        frames = psi.shape[0], interval = 10, blit = False)

    if (filename_save is not None):
        mpeg_writer = animation.FFMpegWriter(fps = 24, bitrate = 10000,
            codec = "libx264", extra_args = ["-pix_fmt", "yuv420p"])
        anim.save("{}.mp4".format(filename_save), writer = mpeg_writer)

    return anim

def animate_contourf(x, y, t, psi, filename_save = None):
    """Function that generates an animation of psi over time in a colored contour plot."""
    plt.style.use("ggplot")
    fig = plt.figure(figsize = (8, 6))
    ax = fig.add_subplot(1, 1, 1)
    
    plot_title = "Streamfunction $\psi(x,y,t)$ at time {}"
    ax.set_title(plot_title.format(t[0]), fontname = "serif", fontsize = 17)
    ax.set_xlabel("x [dim-less]", fontname = "serif", fontsize = 12)
    ax.set_ylabel("y [dim-less]", fontname = "serif", fontsize = 12)
    print(psi.shape)
    CS = ax.contourf(x, y, psi[0, :, :], 20, vmin = -1, vmax = 1, cmap = plt.cm.RdBu_r)
    fig.colorbar(CS, orientation = "vertical")

    def update_frame(frame):
        """Function for updating frame in animation."""
        ax.clear()
        ax.set_title(plot_title.format(t[frame]),
            fontname = "serif", fontsize = 16)
        ax.contourf(x, y, psi[frame, :, :], 20, vmin = -1, vmax = 1, cmap = plt.cm.RdBu_r)

        return ax

    anim = animation.FuncAnimation(fig, update_frame,
        frames = psi.shape[0], interval = 10, blit = False)

    if (filename_save is not None):
        mpeg_writer = animation.FFMpegWriter(fps = 24, bitrate = 10000,
            codec = "libx264", extra_args = ["-pix_fmt", "yuv420p"])
        anim.save("{}.mp4".format(filename_save), writer = mpeg_writer)

    return anim


filename_2D = "../build-Project5/results_2d.txt"
x, y, t, psi = get_data(filename_2D)

psi_0 = psi.reshape((psi.shape[0], 41, 41))

#contour_results(x, y, psi_0)
#anim = animate_pcolormesh(x, y, t, psi_0)
anim = animate_contourf(x, y, t, psi_0)

plt.show()