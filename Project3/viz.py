
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

#with open("../build-Project3/results.txt", "r") as file_result:
#	for line_num, linestring in enumerate(file_result):

filename = "../build-Project3/{}.txt"

data_earth = np.loadtxt(filename.format("earth"), dtype = np.float64, skiprows = 1)
data_jupiter = np.loadtxt(filename.format("jupiter"), dtype = np.float64, skiprows = 1)
data_mars = np.loadtxt(filename.format("mars"), dtype = np.float64, skiprows = 1)
print(data_earth.shape)

# plt.figure(figsize = (8, 8))
# plt.axis("equal")
# plt.plot(data_earth[:, 1], data_earth[:, 2])
# plt.plot(data_jupiter[:, 1], data_jupiter[:, 2])
# plt.plot(data_mars[:, 1], data_mars[:, 2])
# plt.show()

fig = plt.figure(figsize = (10, 7))
fig.set_size_inches(8, 8)

ax = plt.axes(xlim=(-8, 8), ylim=(-8, 8))
ax.set_aspect("equal")

sun = plt.Circle((0,0, 0.0), 0.5, fc = "y")
earth = plt.Circle((data_earth[0, 1], data_earth[0, 2]), 0.04, fc='g')
mars = plt.Circle((data_mars[0, 1], data_mars[0, 2]), 0.03, fc='r')
jupiter = plt.Circle((data_jupiter[0, 1], data_jupiter[0, 2]), 0.2, fc='orange')

def init():
    earth.center = (data_earth[0, 1], data_earth[0, 2])
    jupiter.center = (data_jupiter[0, 1], data_jupiter[0, 2])
    mars.center = (data_mars[0, 1], data_mars[0, 2])
    ax.add_patch(sun)
    ax.add_patch(earth)
    ax.add_patch(mars)
    ax.add_patch(jupiter)
    return []

def animation_manage(i,earth,mars,jupiter):
    animate_earth(i,earth)
    animate_mars(i,mars)
    animate_jupiter(i,jupiter)
    return []

def animate_earth(i, patch):
    x = data_earth[i, 1]
    y = data_earth[i, 2]
    patch.center = (x, y)
    return patch,

def animate_mars(i, patch):
    x = data_mars[i, 1]
    y = data_mars[i, 2]
    patch.center = (x, y)
    return patch,

def animate_jupiter(i, patch):
    x = data_jupiter[i, 1]
    y = data_jupiter[i, 2]
    patch.center = (x, y)
    return patch,


anim = animation.FuncAnimation(fig, animation_manage, 
                               init_func=init, 
                               fargs=(earth,mars,jupiter,),
                               frames=data_earth.shape[0], 
                               interval=10,
                               blit=True)
plt.show()




