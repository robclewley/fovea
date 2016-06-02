from graphics import *
from numpy import *
from matplotlib.pyplot import *
import time

x = linspace(-5, 5)
y = x ** 2

plotter = Plotter()
gui = diagnosticGUI(plotter)

plotter.add_fig('Master', title='x**2', xlabel='x', ylabel='y', domain=[(-100, 100), (-100, 100)])
plotter.add_fig('Minor', title='x**3', xlabel='x', ylabel='y', domain=[(-100, 100), (-100, 100)])

plotter.add_layer('x2_data', figure='Master')

plotter.add_data([x,y], layer='x2_data', figure='Master', style='g-')

plotter.add_layer('x3_data', figure='Minor')

x = linspace(-5, 5)
y = x ** 3

plotter.add_data([x, y], figure='Minor', layer='x3_data', style='r-')
plotter.arrange_fig([1,1], {'11':
                                {'name': 'Plot of x**2',
                                 'scale': [(-20,20),(-50,50)],
                                 'layers': ['x2_data'],
                                 'axes_vars': ['x', 'y'],
                                }}, figure='Master')

plotter.arrange_fig([1,1], {'11':
                                {'name': 'Plot of x**3',
                                 'scale': [(-20,20),(-50,50)],
                                 'layers': ['x3_data'],
                                 'axes_vars': ['x', 'y'],
                                }}, figure='Minor')

plotter.auto_scale_domain(figure='Master', subplot='11')
plotter.auto_scale_domain(figure='Minor', subplot='11')

gui.build_plotter((20, 20), with_times=False)

time.sleep(10)

plotter.clear_fig_data('Minor')

gui.build_plotter((20, 20), with_times=False)
