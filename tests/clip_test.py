import PyDSTool as dst
from PyDSTool.Toolbox.phaseplane import Point2D
from matplotlib import pyplot as plt
import numpy as np

from fovea.graphics import *

plotter.domain={'x': [-1,1], 'y': [-2,2]}
plotter.coords=('x','y')

a, b = force_line_to_extent(Point2D((0.25,0)), Point2D((0.5,0.1)),
                            plotter.domain, plotter.coords)
plt.plot(np.array((a, b)).T[0], np.array((a, b)).T[1], 'g')

cc, dd = Point2D((3,-2.5)), Point2D((-5.1,1.8))
plt.plot(np.array((cc, dd)).T[0], np.array((cc, dd)).T[1], 'b:')
c, d = force_line_to_extent(Point2D((3,-2.5)), Point2D((-5.1,1.8)),
                            plotter.domain, plotter.coords)
plt.plot(np.array((c, d)).T[0], np.array((c, d)).T[1], 'r')

e, f = force_line_to_extent(np.array((0.25,0)), np.array((1.5,0.1)),
                            plotter.domain, plotter.coords)
plt.plot(np.array((e, f)).T[0], np.array((e, f)).T[1], 'g')


plt.show()
