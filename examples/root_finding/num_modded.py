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
import numpy

global dm
dm = diagnostic_manager('rootfinding')
plotter = gui.plotter
# -----------------------------------------

def plot_pt(x, f, n, name, col, layer_root, marker='o'):
    """
    Internal diagnostic helper function.
    """
    new_layer = layer_root+'_data_%d'%n
    pt = pp.Point2D(x, f(x))
    plotter.add_point(pt, style=col+marker, name=name+'_%d'%n,
                     layer=new_layer, log=dm.log)
    fs = plotter.figs['Master']
    ax = fs.arrange['11']['axes_obj']
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_ext = xlim[1]-xlim[0]
    y_ext = ylim[1]-ylim[0]
    # add marker text at a distance proportional to size of
    # current axis extent
    plotter.add_text(pt[0]-0.05*x_ext, pt[1]+0.02*y_ext, name,
                    use_axis_coords=False, name=name+'_%d'%n,
                    layer=layer_root+'_text_%d'%n, style=col)

    return pt


def bisection(f, a, b, TOL=0.001, NMAX=100):
    """
    Takes a function f, start values [a,b], tolerance value(optional) TOL and
    max number of iterations(optional) NMAX and returns the root of the equation
    using the bisection method.
    """
    n=1
    dm.log.msg('Call args', a=a, b=b, TOL=TOL, NMAX=NMAX)
    plotter.add_text(0.1, 0.95, 'n=%d'%n, use_axis_coords=True,
                    name='n_value', layer='meta_data', style='k')
    while n<=NMAX:
        dm.log = dm.log.bind(n=n)
        plotter.set_text('n_value', 'n=%d'%n, 'meta_data')
        if n == 1:
            rebuild = True
        else:
            plotter.toggle_display(layer='bisect_text_%d'%(n-1))
            plotter.toggle_display(layer='bisect_data_%d'%(n-1))
            rebuild = False
        ##ISSUE: Must include subplot arg to ensure layer ends up in fig_struct.arrange. This is unexpected and inconvenient.
        plotter.add_layer('bisect_data_%d'%n, subplot= '11')
        plotter.add_layer('bisect_text_%d'%n, kind='text', subplot= '11')
        c = (a+b)/2.0
        dm.log.msg('Bisect loop', a=a, b=b, c=c)
        a_pt = plot_pt(a, f, n, 'a', 'r', 'bisect', 'o')
        b_pt = plot_pt(b, f, n, 'b', 'g', 'bisect', 'o')
        c_pt = plot_pt(c, f, n, 'c', 'k', 'bisect', 'x')

        if f(c)==0 or (b-a)/2.0 < TOL:
            dm.log.msg('Success', fval=f(c), err=(b-a)/2.0)
            dm.log = dm.log.unbind('n')
            plotter.show(rebuild=rebuild)
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
            plotter.show(rebuild=rebuild)
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
    dm.log.msg('Call args', x0=x0, x1=x1, TOL=TOL, NMAX=NMAX)
    plotter.add_text(0.1, 0.95, 'n=%d'%n, use_axis_coords=True,
                    name='n_value', layer='meta_data', style='k')

    while n<=NMAX:
        dm.log = dm.log.bind(n=n)
        plotter.set_text('n_value', 'n=%d'%n, 'meta_data')
        if n == 1:
            rebuild = True
        else:
            plotter.toggle_display(layer='secant_text_%d'%(n-1))
            plotter.toggle_display(layer='secant_data_%d'%(n-1))
            rebuild = False
        plotter.add_layer('secant_data_%d'%n, subplot= '11')
        plotter.add_layer('secant_text_%d'%n, kind='text', subplot= '11')
        x2 = x1 - f(x1)*((x1-x0)/(f(x1)-f(x0)))
        x0_pt = plot_pt(x0, f, n, 'x0', 'r', 'secant', 'o')
        x1_pt = plot_pt(x1, f, n, 'x1', 'g', 'secant', 'o')
        x2_pt = plot_pt(x2, f, n, 'x2', 'k', 'secant', 'x')
        fs = plotter.figs['Master']
        fx2 = f(x2)
        xleft = fs.domain[0][0]
        xright = fs.domain[0][1]
        gradient = (x1_pt[1]-x0_pt[1])/(x1_pt[0]-x0_pt[0])
        yleft = gradient*(xleft-x1_pt[0]) + x1_pt[1]
        yright = gradient*(xright-x1_pt[0]) + x1_pt[1]
        plotter.add_line_by_points((pp.Point2D(xleft, yleft),
                                 pp.Point2D(xright, yright)),
                                layer='secant_data_%d'%n,
                                name='secantline_%d'%n, style='b-', log=dm.log)
        plotter.add_line_by_points((pp.Point2D(x2, 0),
                                 pp.Point2D(x2, fx2)),
                                layer='secant_data_%d'%n,
                                name='vertline_%d'%n, style='r-', log=dm.log)
        dm.log.msg('Secant loop', x0=x0, x1=x1, x2=x2)
        if abs(x2-x1) < TOL:
            dm.log.msg('Success', fval=fx2, err=abs(x2-x1))
            dm.log = dm.log.unbind('n')
            plotter.show(rebuild=rebuild)
            return x2
        else:
            # The missing increment was discovered when
            # layer with same name ('secant_data_1') tried to be
            # recreated!
            n = n+1
            dm.log.msg('Step: x1->x0, x2->x1', err=abs(x2-x1))
            x0 = x1
            x1 = x2
            plotter.show(rebuild=rebuild)
    dm.log.msg('Failure', status='fail', fval=f(x1), err=abs(x2-x1))
    dm.log = dm.log.unbind('n')
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
