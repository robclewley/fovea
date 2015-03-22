"""
polygon testing
"""
from __future__ import division
import PyDSTool.Toolbox.phaseplane as pp
import PyDSTool as dst

import matplotlib as mpl
import matplotlib.pyplot as plt
import shapely.geometry as geom
from shapely.geometry import polygon as P
from descartes.patch import PolygonPatch

def plot_coords(ax, ob, col, style, lw=1):
    x, y = ob.xy
    ax.plot(x, y, style+col, lw=lw, zorder=1)
    plt.draw()

def test_line1(x):
    y = 3*x+2*math.pi
    return y

def test_func1(pt):
    x, y = pt
    return y - (test_line1(x))


def test_line2(x):
    y = 1.2*x-3.5*math.pi
    return y

def test_func2(pt):
    x, y = pt
    return y - (test_line2(x))

# TEST 1
def test1():
    poly1 = P.Polygon([(1,2), (0,0), (0.5,-2), (2.5,-1), (2,1)])
    poly2 = P.Polygon([(1,1), (2,-3), (6,-1), (5,1), (3,3)])

    cent1 = poly1.centroid
    cent2 = poly2.centroid

    assert not poly1.within(poly2)
    assert not poly2.within(poly1)

    poly3 = poly2.union(poly1)
    inter = poly2.intersection(poly1)

    plt.figure()
    ax = plt.gca()
    plot_coords(ax, poly1.exterior, 'r', 'o-', lw=3)
    plot_coords(ax, poly2.exterior, 'b', 'o-', lw=3)
    ax.plot(cent1.x, cent1.y, 'go')
    ax.plot(cent2.x, cent2.y, 'gs')
    plot_coords(ax, poly3.exterior, 'k', 'o-')

    patch = PolygonPatch(inter, facecolor='#dd3333', alpha=0.5, zorder=0.5)
    ax.add_patch(patch)
    plt.show()

# TEST 2

import numpy as np
import math
import domain2D as dom

xs = np.linspace(-6,12,2)
ys1 = test_line1(xs)
ys2 = test_line2(xs)

plt.figure()
ax = plt.gca()
ax.set_aspect('equal')

plt.plot(xs,ys1, 'k-', lw=4)
plt.plot(xs,ys2, 'k-', lw=4)

obj = dom.polygon_domain(pp.Point2D(1,0.5), pp.Point2D(1,-1),
                     [test_func1, test_func2], max_radius_fac=20,
                     edge_len_rtol=4)

plot_coords(ax, obj.polygon.exterior, '', 'o-', lw=3)
plt.show()

for i in range(60):
    plt.draw()
    obj.grow_step()
    plot_coords(ax, obj.polygon.exterior, '', '-')

plot_coords(ax, obj.polygon.exterior, 'k', 'o')
