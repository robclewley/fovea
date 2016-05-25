from fovea.graphics import *
from numpy import *
from matplotlib.pyplot import *

x = linspace(-5, 5)
y = x**2

plotter = Plotter()
gui = diagnosticGUI(plotter)

plotter.add_fig('Master', title='quadratic', xlabel='x', ylabel='y', domain=[(-5, 5), (0,25)])

plotter.add_layer('fn_data')

plotter.add_data([x,y], layer='fn_data', style='g-')

plotter.arrange_fig([1,1], {'11':
                           {'name': 'Plot of x**2',
                            'scale': [(-5,5),(0,25)],
                            'layers': ['fn_data'],
                            'axes_vars': ['x', 'y'],
                           }})

gui.build_plotter((14, 6), with_times=False)
