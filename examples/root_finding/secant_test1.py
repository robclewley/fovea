import math

from num_modded import secant, plotter, gui, dm
from fovea import *


from test_funcs import f1


# Boilerplate

# domain of interest
DOI = ([0,2.1],[-2.5,2])

plotter.clean()
plotter.addFig('Master',
               title='Root Finding Diagnostic',
               xlabel='x', ylabel='y',
               domain=DOI)

plotter.addLayer('fn_data')
plotter.addLayer('meta_data', kind='text')

plotter.arrangeFig([1,1], {'11':
                           {'name': 'Secant method',
                            'scale': DOI,
                            'layers': '*',  # all layers will be selected
                            'axes_vars': ['x', 'y']}
                           })

gui.buildPlotter2D((8,8), with_times=False)

plotter.addHLine(0, style='k:', layer='fn_data')
plotter.addVLine(0, style='k:', layer='fn_data')
xs = npy.linspace(DOI[0][0], DOI[0][1], 500)
ys = f1(xs)
plotter.addData([xs, ys], layer='fn_data', style='k-')
plt.show()

root = secant(f1, 1.05, 2)

# demonstration of database log
db = dm.log.get_DB()
dblist = db.all()
n_max = dblist[-1]['n']
