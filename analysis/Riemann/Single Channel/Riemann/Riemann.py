import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pylab
import numpy as np
import mpmath
mpmath.dps = 5

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
import numpy as np

if sys.argv[1] == 'both':
  real = True
  imag = True
  col = 'r'
  axlabel = "Domain"
elif sys.argv[1] == 'real':
  real = True
  imag = False
  col = 'b'
  axlabel = "Range real"
elif sys.argv[1] == 'imag':
  real = False
  imag = True
  col = 'b'
  axlabel = "Range imag"
  
if sys.argv[2] == 'both':
  positive = True
  negative = True
elif sys.argv[2] == 'pos':
  positive = True
  negative = False
elif sys.argv[2] == 'neg':
  positive = False
  negative = True

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-5, 5, 0.05)
Y = np.arange(-5, 5, 0.05)
X, Y = np.meshgrid(X, Y)
  
fr = lambda x:x.real
fi = lambda x:x.imag
if positive:
  if real:
    Zpos = map(fr, np.sqrt(X+Y*1.0j))
    ax.plot_wireframe(X, Y, Zpos, rstride=10, cstride=10)
  if imag:
    Zpos = map(fi, np.sqrt(X+Y*1.0j))
    ax.plot_wireframe(X, Y, Zpos, rstride=10, cstride=10, color=col)
if negative:
  if real:
    Zneg = map(fr, -np.sqrt(X+Y*1.0j))
    ax.plot_wireframe(X, Y, Zneg, rstride=10, cstride=10)
  if imag:
    Zneg = map(fi, -np.sqrt(X+Y*1.0j))
    ax.plot_wireframe(X, Y, Zneg, rstride=10, cstride=10, color=col)
  
'''
theta = np.linspace(0, 4*np.pi, 100)
z = np.exp((theta/2.0+np.pi)*1.0j)
x = np.cos(theta)
y = np.sin(theta)
ax.plot(x, y, z, label='parametric curve', c='r')

vals = np.linspace(-6.0, 0.0, 100)
z = map(fi, -np.sqrt(vals))
x = vals
y = np.linspace(0.0, 0.0, 100)
ax.plot(x, y, z, label='parametric curve', c='r')
'''

ax.set_xlabel('Domain real')
ax.set_ylabel('Domain imag')
ax.set_zlabel(axlabel)  
pylab.show()