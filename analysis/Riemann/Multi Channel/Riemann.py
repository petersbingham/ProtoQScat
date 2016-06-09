import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pylab
import numpy as np
import sympy.mpmath as mpmath
mpmath.dps = 5

REMOVEDISCONT = True

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
import numpy as np

if sys.argv[1] == 'both':
  real = True
  imag = True
  axlabel = "p"
elif sys.argv[1] == 'real':
  real = True
  imag = False
  axlabel = "p real"
elif sys.argv[1] == 'imag':
  real = False
  imag = True
  axlabel = "p imag"
  
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

XX = np.arange(-5, 5, 0.05)
if REMOVEDISCONT:
  Ys = [np.arange(-5, 0, 0.05), np.arange(0, 5, 0.05)]
else:
  Ys = [np.arange(-5, 5, 0.05)]
  
fr = lambda x:x.real
fi = lambda x:x.imag

def _val(X,Y,Th):
  return np.sqrt((X+Y*1.0j)-Th)

desc = {0.0:'r',3.0:'b'}  
for Th in sorted(desc.keys()):
  for YY in Ys:
    X, Y = np.meshgrid(XX, YY)
    if positive:
      if real:
        Zpos = map(fr, _val(X,Y,Th))
        ax.plot_wireframe(X, Y, Zpos, rstride=10, cstride=10, color=desc[Th])
      if imag:
        Zpos = map(fi, _val(X,Y,Th))
        ax.plot_wireframe(X, Y, Zpos, rstride=10, cstride=10, color=desc[Th])
    if negative:
      if real:
        Zneg = map(fr, -_val(X,Y,Th))
        ax.plot_wireframe(X, Y, Zneg, rstride=10, cstride=10, color=desc[Th])
      if imag:
        Zneg = map(fi, -_val(X,Y,Th))
        ax.plot_wireframe(X, Y, Zneg, rstride=10, cstride=10, color=desc[Th])
  
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

ax.set_xlabel('E real')
ax.set_ylabel('E imag')
ax.set_zlabel(axlabel)  
pylab.show()