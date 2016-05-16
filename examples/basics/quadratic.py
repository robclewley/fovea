from fovea.graphics import *
from numpy import *
from matplotlib.pyplot import *

x = linspace(-5, 5)
y = x**2

plotter = plotter2D()
gui = diagnosticGUI(plotter)

plotter.addFig('Master', title='quadratic', xlabel='x', ylabel='y', domain=[(-5, 5), (0,25)])

plotter.addLayer('fn_data')

plotter.addData([x,y], layer='fn_data', style='g-')

plotter.arrangeFig([1,1], {'11':
                           {'name': 'Plot of x**2',
                            'scale': [(-5,5),(0,25)],
                            'layers': ['fn_data'],
                            'axes_vars': ['x', 'y'],
                           }})

gui.buildPlotter2D((14, 6), with_times=False)
