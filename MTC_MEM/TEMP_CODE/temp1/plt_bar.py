import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

path = r'C:\Users\13323\Downloads\res.txt'
data = np.loadtxt(path)
data = data.T
X,Y = np.meshgrid(data[0, :],data[1,:])
plt.figure(dpi=100, figsize=(8,6))
plt.title('Steady state solution for CO$_2$\n'
         "30 μm electrode on left and 60 μm electrode on right")
plt.contour(X,Y,data[2,:],cmap='plasma')
#c.set_cmap("plasma")
#plot(mesh, linewidth=0.1, zorder=1)
plt.colorbar(fraction=0.015, pad=0.08)
plt.gca().ticklabel_format(style='sci', scilimits=(-1,1))
plt.xlabel('distance (m)')

plt.show()