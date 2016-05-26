import math

from num_modded import secant, plotter, gui, dm
from fovea import *

from test_funcs import f2

# domain of interest
DOI = ([0,20],[-1,0.5])

plotter.clean()
plotter.add_fig('Master',
               title='Root Finding Diagnostic',
               xlabel='x', ylabel='y',
               domain=DOI)

dm.use_dir('secant2_success')
dm.make_log('log.json')
plotter.dm = dm
# default to wait for user input on each iteration
plotter.wait_status = True

plotter.add_layer('fn_data')
plotter.add_layer('meta_data', kind='text')

plotter.arrange_fig([1,1], {'11':
                           {'name': 'Secant method',
                            'scale': DOI,
                            # "*" functionality is broken
                            #'layers': '*',  # all layers will be selected
                            'layers': ['fn_data', 'meta_data'],
                            'axes_vars': ['x', 'y']}
                           })

gui.build_plotter((8,8), with_times=False)

plotter.add_hline(0, style='k:', layer='fn_data')
plotter.add_vline(0, style='k:', layer='fn_data')
xs = npy.linspace(DOI[0][0], DOI[0][1], 1000)
ys = f2(xs)
plotter.add_data([xs, ys], layer='fn_data', style='k-')

root = secant(f2, 3, 12)
print("Root is x=%.4f"%root)

# demonstration of database log
db = dm.log.get_DB()
dblist = db.all()
n_max = dblist[-1]['n']
