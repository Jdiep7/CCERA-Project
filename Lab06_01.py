import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, sin

nRows, nCols = 500, 500
mapData = 0.3*np.ones((nRows,nCols))   # gray background

zMax = 1.
k = 2.
x = np.linspace(-1.,1.,nRows)
y = np.linspace(-1.,1.,nCols)
for i in range(nCols) :
    for j in range(nRows) :
        r = sqrt(x[j]*x[j] + y[i]*y[i])
        mapData[j][i] = sin(k*r)/ r 

fig = plt.figure(figsize=(10.5,8.))
ax = fig.add_subplot(111)
ax.set_title("sin(kr)/r")
ax.set_xlabel("x")
ax.set_ylabel("y")

ax.patch.set_facecolor('white')

#show the plot
im = ax.imshow(mapData,extent=[-10,10.,-10,10.],aspect=1.0)
plt.colorbar(im, use_gridspec=True)
plt.show()

    