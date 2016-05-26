from graphics import *
from numpy import *
from matplotlib.pyplot import *

x = linspace(-5, 5)
y = x ** 2

plotter = Plotter()
gui = diagnosticGUI(plotter)

plotter.add_fig('Master', title='domain_test', xlabel='x', ylabel='y', domain=[(-100, 100), (-100, 100)])

plotter.add_layer('fn_data')

plotter.add_data([x,y], layer='fn_data', style='g-')


plotter.add_layer('layer_dos')
x = linspace(-5, 5)
y = x ** 3
plotter.add_data([x, y], layer='layer_dos', style='r-')
plotter.arrange_fig([1,1], {'11':
                           {'name': 'Plot of x**2',
                            'scale': [(-20,20),(-50,50)],
                            'layers': ['fn_data', 'layer_dos'],
                            'axes_vars': ['x', 'y'],
                           }})

plotter.auto_scale_domain(subplot='11')

print(plotter.active_layer)
# ('Master', 'layer_dos')
plotter.set_active_layer('fn_data')
print(plotter.active_layer)
# ('Master', 'fn_data')

gui.build_plotter((20, 20), with_times=False)
