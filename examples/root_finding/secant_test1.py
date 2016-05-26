import math

from num_modded import secant, plotter, gui, dm
from fovea import *


from test_funcs import f1

# domain of interest
DOI = ([0,2.1],[-2.5,2])
#DOI = ([0,6],[-2.5,30])

plotter.clean()
plotter.add_fig('Master',
               title='Root Finding Diagnostic',
               xlabel='x', ylabel='y',
               domain=DOI)

dm.use_dir('secant1_fail2')
dm.make_log('log.json')
plotter.dm = dm
# default to wait for user input on each iteration
plotter.wait_status = True

plotter.add_layer('fn_data')
plotter.add_layer('meta_data', kind='text')

plotter.arrange_fig([1,1], {'11':
                           {'name': 'Secant method',
                            'scale': DOI,
                            'layers': '*',  # all layers will be selected
                            'axes_vars': ['x', 'y']}
                           })

gui.build_plotter((8,8), with_times=False)

plotter.add_hline(0, style='k:', layer='fn_data')
plotter.add_vline(0, style='k:', layer='fn_data')
xs = npy.linspace(DOI[0][0], DOI[0][1], 500)
ys = f1(xs)
plotter.add_data([xs, ys], layer='fn_data', style='k-')

root = secant(f1, 1.05, 2)
print("Root is x=%.4f"%root)

# demonstration of database log
db = dm.log.get_DB()
dblist = db.all()
n_max = dblist[-1]['n']
