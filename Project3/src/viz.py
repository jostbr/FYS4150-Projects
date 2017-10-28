
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

# ========================= GLOBAL VARIABLES =========================
file_directory = "../build-Project3/"   # Path to result data files

#planet_masses = {"sun": 2.0E+30, "earth": 6.0E+24, "jupiter": 1.9E+27, "mars": 6.6E+23,
#    "venus": 4.9E+24, "saturn": 5.5E+26, "mercury": 3.3E+23,
#    "uranus": 8.8E+25, "neptune": 1.03E+26, "pluto": 1.31E+22, "moon": 7.34E+22}
circle_sizes = {"sun": 0.1, "earth": 0.08, "jupiter": 0.4, "mars": 0.04, "venus": 0.1,
    "saturn": 0.2, "mercury": 0.02, "uranus": 0.1, "neptune": 0.12, "pluto": 0.14,
    "moon": 0.008}
circle_colors = {"sun": "#d6c200", "earth": "#2666FE", "jupiter": "#f99805", "mars": "#f24519",
    "venus": "#c56300", "saturn": "#DAA520", "mercury": "#965402", "uranus": "#A6D5CA",
    "neptune": "#4D8FAC", "pluto": "#8b8b8b", "moon": "#000000"}

planets = list(circle_sizes.keys())     # Names of all planets
data = dict()   # To hold data for all planets

for planet in planets:
    filename = planet + ".txt"              # Filename of data file
    file_path = file_directory + filename   # Corresponding path

    if (os.path.exists(file_path)):
        data[planet] = np.loadtxt(file_path, dtype = np.float64, skiprows = 1)  # Load in data from file

some_planet = list(data.keys())[0]  # When needing to access an arbitrary element in data
# ====================================================================

# Function top plot trajectories of all computed planets.
def plot_trajectories_2D(data):
    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize = (8, 8))
    ax.set_aspect("equal")

    if ("sun" not in list(data.keys())):
        sun = plt.Circle((0.0, 0.0), 0.03, fc = "#d6c200")    # Need some sun
        ax.add_patch(sun)

    for planet in data.keys():
        ax.plot(data[planet][:, 1], data[planet][:, 2], color = circle_colors[planet],
            label = planet.title())

    #ax.text(0.03, 0.03, "Sun")
    ax.set_title("Planet trajectories over {:.2f} years".format(data[some_planet][-1, 0]),
        fontname = "serif", fontsize = 20)
    ax.set_xlabel("x [AU]", fontname = "serif", fontsize = 12)
    ax.set_ylabel("y [AU]", fontname = "serif", fontsize = 12)
    ax.legend()

def plot_trajectories_3D(data):
    #plt.style.use("default")
    fig = plt.figure(figsize = (11, 7))
    ax = Axes3D(fig)

    end_time_index = -1

    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
    x = 0.05*np.cos(u)*np.sin(v)
    y = 0.05*np.sin(u)*np.sin(v)
    z = 0.05*np.cos(v)

    # Plot the surface
    #ax.plot_surface(x, y, z, color = "y")
    #ax.text(0.05, 0.05, 0.05, "Sun")

    for planet in data.keys():
        ax.plot(data[planet][:end_time_index, 1], data[planet][:end_time_index, 2],
            data[planet][:end_time_index, 3], color = circle_colors[planet], label = planet.title())

    #ax.set_zlim(-1, 1)
    #ax.set_xlim(-1,1)
    #ax.set_ylim(-1,1)


    ax.set_title("Planet trajectories over {:.2f} years".format(data[some_planet][end_time_index, 0]),
        fontname = "serif", fontsize = 18)
    ax.set_xlabel("x [AU]", fontname = "serif", fontsize = 16)
    ax.set_ylabel("y [AU]", fontname = "serif", fontsize = 16)
    ax.set_zlabel("z [AU]", fontname = "serif", fontsize = 16)
    ax.legend()

# Function to animate trajectories of all computed planets
def animate_trajectories(data):
    plt.style.use("default")
    fig = plt.figure(figsize = (10, 7))
    fig.set_size_inches(8, 8)

    ax = plt.axes(xlim=(-20, 20), ylim=(-20, 20))
    ax.set_aspect("equal")
    ax.set_title("Planet positions after {:.2f} years".format(data[some_planet][0, 0]),
        fontname = "serif", fontsize = 18)
    ax.set_xlabel("x-coordinate [AU]", fontname = "serif", fontsize = 12)
    ax.set_ylabel("y-coordinate [AU]", fontname = "serif", fontsize = 12)

    circles = dict()    # To hold plt.Circle objects representing planets
    lines = list()      # Only be able to create legend with circle markers
    texts = dict()

    for planet in data.keys():
        texts[planet] = ax.text(data[planet][0, 1], data[planet][0, 2], planet.title(), fontsize = 6)

    if ("sun" not in list(data.keys())):
        sun = plt.Circle((0.0, 0.0), 0.1, fc = "#d6c200")    # Need some sun
        ax.add_patch(sun)

    for planet in data.keys():
        circles[planet] = plt.Circle((data[planet][0, 1], data[planet][0, 2]),
            circle_sizes[planet], fc = circle_colors[planet])   # Initiate circles
        lines.append(plt.Line2D(range(1), range(1), color = "none", marker = "o",
            markerfacecolor = circle_colors[planet], label = planet.title()))   # Only for legend below

    # Initialize animation
    def init():
        for planet in circles.keys():
            circles[planet].center = (data[planet][0, 1], data[planet][0, 2])
            ax.add_patch(circles[planet])

        return []

    # Make sure each planet/circle position is updated
    def animation_manage(i, planets):
        for planet in circles.keys():
            animate_planet(i, planet)
            animate_text(i, planet)

        ax.set_title("Planet positions after {:.2f} years".format(data[some_planet][i, 0]),
            fontname = "serif", fontsize = 18)
        return []

    # Update each planet/circle position.
    def animate_planet(i, planet):
        x = data[planet][i, 1]
        y = data[planet][i, 2]
        circles[planet].center = (x, y)
        return circles[planet]

    def animate_text(i, planet):
        texts[planet].set_x(data[planet][i, 1] + 0.1)
        texts[planet].set_y(data[planet][i, 2] + 0.1)
        return texts[planet],

    fig.tight_layout()
    ax.legend(lines, [planet.title() for planet in data.keys()])     # Label all planets
    anim = animation.FuncAnimation(fig, animation_manage, init_func = init,
        fargs=(list(circles.keys()),), frames = data[some_planet].shape[0],
        interval = 10, blit = True)

    return anim

plot_trajectories_2D(data)
#plot_trajectories_3D(data)
#anim = animate_trajectories(data)

plt.show()




