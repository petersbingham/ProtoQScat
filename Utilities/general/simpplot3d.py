import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import pylab
import numpy as np
import sympy.mpmath as mpmath
mpmath.dps = 5


def plot(Rs, Is, fun, real, indepLabel="", depLabel=""):
    real = True
    
    if real:
        axlabel = depLabel + " real"
    else:
        axlabel = depLabel + " imag"
      
    fig = plt.figure()
    ax = fig.gca(projection='3d')
      
    fr = lambda x:x.real
    fi = lambda x:x.imag

    Rs, Is = np.meshgrid(Rs,Is)
    if real:
        Zs = map(fr, fun(Rs,Is))
    else:
        Zs = map(fi, fun(Rs,Is))
    ax.plot_wireframe(Rs, Is, Zs, rstride=10, cstride=10)
    
    ax.set_xlabel(indepLabel+' real')
    ax.set_ylabel(indepLabel+' imag')
    ax.set_zlabel(axlabel)  
    pylab.show()


if __name__ == '__main__':
    Rs = np.arange(-5, 5, 0.05)
    Is = np.arange(-5, 5, 0.05)

    s = plot(Rs, Is, lambda R,I:np.sqrt(R+I*1.0j), True)