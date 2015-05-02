import math

from num_modded import bisection, plotter, gui, dm
from fovea import *

def f1(x):
    """
    Test function 1
    """
    return x*x*x - math.pi*x + math.e/100


# Boilerplate

# domain of interest
DOI = ([-1,1],[-2,2])

plotter.clean()
plotter.addFig('Master',
               title='Root Finding Diagnostic',
               xlabel='x', ylabel='y',
               domain=DOI)

plotter.addLayer('bisect_data')


plotter.arrangeFig([1,1], {'11': {'name': 'plane',
                                  'scale': DOI,
                                  'layers': ['bisect_data'],
                                  'axes_vars': ['x', 'y']}
                           })

# Redundant: last added layer becomes active by default
#plotter.set_active_layer('bisect_data')

gui.buildPlotter2D((8,8), with_times=False)

root = bisection(f1, 0.1, 1)