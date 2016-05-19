import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pylab
import numpy as np
import sympy.mpmath as mpmath
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
  zlabel = "p"
elif sys.argv[1] == 'real':
  real = True
  imag = False
  col = 'b'
  zlabel = "p real"
elif sys.argv[1] == 'imag':
  real = False
  imag = True
  col = 'b'
  zlabel = "p imag"
  
if sys.argv[2] == 'both':
  physical = True
  nonphysical = True
elif sys.argv[2] == 'physical':
  physical = True
  nonphysical = False
elif sys.argv[2] == 'nonphysical':
  physical = False
  nonphysical = True

fig = plt.figure()
ax = fig.gca(projection='3d')
  
fr = lambda x:x.real
fi = lambda x:x.imag

def _plotPhysical(preal):
  X = np.arange(-5, 5, 0.05)
  Y = np.arange(0, 5, 0.05)
  X, Y = np.meshgrid(X, Y)
  Zpos = map(fr if preal else fi, np.sqrt(X+Y*1.0j))
  ax.plot_wireframe(X, Y, Zpos, rstride=10, cstride=10)
  
  X = np.arange(-5, 5, 0.05)
  Y = np.arange(-5, 0, 0.05)
  X, Y = np.meshgrid(X, Y)
  Zneg = map(fr if preal else fi, -np.sqrt(X+Y*1.0j))
  ax.plot_wireframe(X, Y, Zneg, rstride=10, cstride=10)

def _plotNonPhysical(preal):
  X = np.arange(-5, 5, 0.05)
  Y = np.arange(-5, 0, 0.05)
  X, Y = np.meshgrid(X, Y)
  Zpos = map(fr if preal else fi, np.sqrt(X+Y*1.0j))
  ax.plot_wireframe(X, Y, Zpos, rstride=10, cstride=10)
  
  X = np.arange(-5, 5, 0.05)
  Y = np.arange(0, 5, 0.05)
  X, Y = np.meshgrid(X, Y)
  Zneg = map(fr if preal else fi, -np.sqrt(X+Y*1.0j))
  ax.plot_wireframe(X, Y, Zneg, rstride=10, cstride=10)  
    
if real:
  if physical:
    _plotPhysical(True)
  if nonphysical:
    _plotNonPhysical(True)
if imag:
  if physical:
    _plotPhysical(False)
  if nonphysical:
    _plotNonPhysical(False)
  
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
ax.set_zlabel(zlabel) 
pylab.show()