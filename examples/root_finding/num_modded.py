"""
Bisection, Secant & Newton Raphson Method.
**For test purposes only**

Gist at: https://gist.github.com/robclewley/ca3df340d39f6501f124
"""

import math
"""
* Variable Description:
*
* f      : Given function
* f_     : Derivative of f
* [a, b] : End point values
* TOL    : Tolerance
* NMAX   : Max number of iterations
"""

# -----------------------------------------
# Diagnostic imports
import tinydb as tdb
import fovea
import fovea.graphics as gx
from fovea.graphics import tracker
from fovea.graphics import gui
from fovea.diagnostics import diagnostic_manager
from PyDSTool.Toolbox import phaseplane as pp

global dm
dm = diagnostic_manager('rootfinding', 'bisect1_log.json')
plotter = gui.plotter
# -----------------------------------------


def bisection(f, a, b, TOL=0.001, NMAX=100):
    """
    Takes a function f, start values [a,b], tolerance value(optional) TOL and
    max number of iterations(optional) NMAX and returns the root of the equation
    using the bisection method.
    """
    def plot_pt(x, n, name, col, marker='o'):
        """
        Internal diagnostic helper function.
        """
        pt = pp.Point2D(x, f(x))
        plotter.addPoint(pt, style=col+marker, name=name+'_%d'%n,
                         layer='bisect_data_%d'%n, log=dm.log)
        fs = plotter.figs['Master']
        ax = fs.arrange['11']['axes_obj']
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        x_ext = xlim[1]-xlim[0]
        y_ext = ylim[1]-ylim[0]
        # add marker text at a distance proportional to size of
        # current axis extent
        plotter.addText(pt[0]-0.05*x_ext, pt[1]+0.02*y_ext, name,
                        use_axis_coords=False, name=name+'_%d'%n,
                        layer='bisect_text_%d'%n, style=col)
        return pt

    n=1

    plotter.addText(0.1, 0.95, 'n=%d'%n, use_axis_coords=True,
                    name='n_value', layer='meta_data', style='k')
    while n<=NMAX:
        dm.log = dm.log.bind(n=n)
        plotter.setText('n_value', 'n=%d'%n, 'meta_data')
        if n == 1:
            rebuild = True
        else:
            plotter.toggleDisplay(layer='bisect_text_%d'%(n-1))
            plotter.toggleDisplay(layer='bisect_data_%d'%(n-1))
            rebuild = False
        plotter.addLayer('bisect_data_%d'%n)
        plotter.addLayer('bisect_text_%d'%n, kind='text')
        c = (a+b)/2.0
        dm.log.msg('Bisect loop', a=a, b=b, c=c)
        a_pt = plot_pt(a, n, 'a', 'r', 'o')
        b_pt = plot_pt(b, n, 'b', 'g', 'o')
        c_pt = plot_pt(c, n, 'c', 'k', 'x')

        if f(c)==0 or (b-a)/2.0 < TOL:
            dm.log.msg('Success', fval=f(c), err=(b-a)/2.0)
            dm.log = dm.log.unbind('n')
            return c
        else:
            n = n+1
            if f(c)*f(a) > 0:
                dm.log.msg('Same sign')
                a=c
            else:
                dm.log.msg('Opposite sign')
                b=c
            dm.log.msg('Step', err=(b-a)/2.0)
        plotter.show(rebuild=rebuild, wait=True)
    dm.log.msg('Failure', status='fail', fval=f(c), err=(b-a)/2.0)
    dm.log = dm.log.unbind('n')
    return False

def secant(f,x0,x1, TOL=0.001, NMAX=100):
    """
    Takes a function f, start values [x0,x1], tolerance value(optional) TOL and
    max number of iterations(optional) NMAX and returns the root of the equation
    using the secant method.
    """
    n=1
    while n<=NMAX:
        x2 = x1 - f(x1)*((x1-x0)/(f(x0)-f(x1)))
        if x2-x1 < TOL:
            return x2
        else:
            x0 = x1
            x1 = x2
    return False

def newtonraphson(f, f_, x0, TOL=0.001, NMAX=100):
    """
    Takes a function f, its derivative f_, initial value x0, tolerance value(optional) TOL and
    max number of iterations(optional) NMAX and returns the root of the equation
    using the newton-raphson method.
    """
    n=1
    while n<=NMAX:
        x1 = x0 - (f(x0)/f_(x0))
        if x1 - x0 < TOL:
            return x1
        else:
            x0 = x1
    return False

if __name__ == "__main__":

    def func(x):
        """
        Function x^3 - x -2
        We will calculate the root of this function using different methods.
        """
        return math.pow(x,3) - x -2

    def func_(x):
        """
        Derivative of the function f(x) = x^3 - x -2
        This will be used in Newton-Rahson method.
        """
        return 3*math.pow(x,2)-1

    #Invoking Bisection Method
    res = bisection(func,1,2)
    print res

    #Invoking Secant Method
    res = bisection(func,1,2)
    print res

    #Invoking Newton Raphson Method
    res = newtonraphson(func,func_,1)
    print res
