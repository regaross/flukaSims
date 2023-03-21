import mpl_toolkits.mplot3d.axes3d as axes3d
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5


points = h5.File('od_master.h5')['data']['neutron_xyz']
radius = 227.0 # cm

Height = 172 #cm
Radius = 64 # cm

fig = plt.figure(figsize = (10,8))
ax = fig.add_subplot(1, 1, 1, projection='3d')

u = np.linspace(0, np.pi*2, 60)
v = np.linspace(0, np.pi, 100)

x = radius * np.outer(np.cos(u), np.sin(v))
y = radius * np.outer(np.sin(u), np.sin(v))
z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

#Column
ax.plot_surface(x, y, z, color='b', alpha = 0.3)

for point in points:
    ax.scatter(point[0], point[1], point[2], color = 'b')
plt.show()