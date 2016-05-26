import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from fovea.graphics import *

mpl.rcParams['legend.fontsize'] = 10

plotter = Plotter()
gui = diagnosticGUI(plotter)

theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
x = r * np.sin(theta)
y = r * np.cos(theta)

plotter.add_fig('Master', title='3D Curve', xlabel='x', ylabel='y', domain=[(-5, 5), (-5, 5)])
plotter.add_layer('fn_data')
plotter.add_data([x, y, z], layer='fn_data')

plotter.arrange_fig([1,1], {'11':
                           {'name': '3D Spiral Curve',
                            'scale': [(-4, 4),(-4, 5),(-2, 2)],
                            'layers': ['fn_data'],
                            'axes_vars': ['x', 'y', 'z'],
                            'projection': '3d'
                           }})

plotter.auto_scale_domain(subplot='11')

gui.build_plotter((20, 20), with_times=False)
