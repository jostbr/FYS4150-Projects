
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# ========================= GLOBAL VARIABLES =========================
file_directory = "../build-Project3/"   # Path to result data files

planet_masses = {"earth": 6.0E+24, "jupiter": 1.9E+27, "mars": 6.6E+23,
    "venus": 4.9E+24, "saturn": 5.5E+26, "mercury": 3.3E+23,
    "uranus": 8.8E+25, "neptun": 1.03E+26, "pluto": 1.31E+22}
circle_sizes = {"earth": 0.08, "jupiter": 0.4, "mars": 0.04, "venus": 0.1,
    "saturn": 0.2, "mercury": 0.02, "uranus": 0.1, "neptun": 0.12, "pluto": 0.014}
circle_colors = {"earth": "#2666FE", "jupiter": "#f99805", "mars": "#f24519",
    "venus": "#c56300", "saturn": "#DAA520", "mercury": "#965402", "uranus": "#A6D5CA",
    "neptun": "#4D8FAC", "pluto": "#8b8b8b"}

planets = list(circle_sizes.keys())     # Names of all planets
data = dict()   # To hold data for all planets

for planet in planets:
    filename = planet + ".txt"              # Filename of data file
    file_path = file_directory + filename   # Cooresponding path

    if (os.path.exists(file_path)):
        data[planet] = np.loadtxt(file_path,
            dtype = np.float64, skiprows = 1)       # Load in data from file
# ====================================================================

# Function top plot trajectories of all computed planets.
def plot_trajectories(data):
    plt.figure(figsize = (8, 8))
    plt.axis("equal")

    for planet in data.keys():
        plt.plot(data[planet][:, 1], data[planet][:, 2], label = planet.title())

    plt.title("Planet trajectories over {:.2f} years".format(data[list(data.keys())[0]][-1, 0],
        fontname = "serif", fontsize = 14))
    plt.xlabel("x-coordinate [AU]", fontname = "serif", fontsize = 12)
    plt.ylabel("y-coordinate [AU]", fontname = "serif", fontsize = 12)
    plt.legend()

# Function to animate trajectories of all computed planets
def animate_trajectories(data):
    fig = plt.figure(figsize = (10, 7))
    fig.set_size_inches(8, 8)

    ax = plt.axes(xlim=(-8, 8), ylim=(-8, 8))
    ax.set_aspect("equal")

    circles = dict()    # To hold plt.Circle objects representing planets

    sun = plt.Circle((0.0, 0.0), 1.5*max(circle_sizes.values()), fc = "#ffdb00")    # Need some sun
    ax.add_patch(sun)

    for planet in data.keys():
        circles[planet] = plt.Circle((data[planet][0, 1], data[planet][0, 2]),
            circle_sizes[planet], fc = circle_colors[planet])   # Initiate circles

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

        
        return []

    # Update each planet/circle position.
    def animate_planet(i, planet):
        x = data[planet][i, 1]
        y = data[planet][i, 2]
        circles[planet].center = (x, y)
        return circles[planet]


    anim = animation.FuncAnimation(fig, animation_manage, init_func = init,
        fargs=(list(circles.keys()),), frames = data[list(data.keys())[0]].shape[0],
        interval = 10, blit = True)

    return anim

plot_trajectories(data)
anim = animate_trajectories(data)

plt.show()




