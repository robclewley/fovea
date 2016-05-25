import math

from num_modded import bisection, plotter, gui, dm
from fovea import *

from test_funcs import f1

# domain of interest
DOI = ([0,2.1],[-2.5,2])

plotter.clean() # in case rerun in same session
plotter.add_fig('Master',
               title='Root Finding Diagnostic',
               xlabel='x', ylabel='y',
               domain=DOI)

dm.use_dir('bisect1_success')
dm.make_log('log.json')

# connect the diagnostic_manager to the plotter
plotter.dm = dm
# default to wait for user input on each iteration
plotter.wait_status = True

plotter.add_layer('fn_data')
plotter.add_layer('meta_data', kind='text')

plotter.arrange_fig([1,1], {'11':
                           {'name': 'Bisection method',
                            'scale': DOI,
                            #'layers': '*',  # all layers will be selected
                            'layers': ['fn_data', 'meta_data'],
                            'callbacks':'*',
                            'axes_vars': ['x', 'y']}
                           })

gui.build_plotter((8,8), with_times=False)

# permanent axes
plotter.add_hline(0, style='k:', layer='fn_data')
plotter.add_vline(0, style='k:', layer='fn_data')

# plot the function
xs = npy.linspace(DOI[0][0], DOI[0][1], 500)
ys = f1(xs)
plotter.add_data([xs, ys], layer='fn_data', style='k-')

root = bisection(f1, 0.1, 2)
print("Root is x=%.4f"%root)

# demonstration of database log
db = dm.log.get_DB()
dblist = db.all()
n_max = dblist[-1]['n']