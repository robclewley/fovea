from pca_disc import plotter, gui, loopPCA, stretch
from fovea import *

from PyDSTool.Toolbox import synthetic_data as sd

import numpy as np

DOI = [(-10,10),(-10,10)]
dim = 3
d = 2

plotter.clean() # in case rerun in same session
plotter.addFig('Master',
               title='PCA Disc',
               xlabel='x', ylabel='y',
               domain=DOI)

#Plot out several different rotations of the original data.
rot_layers= [['rot_data1','rot_pc1','loD_data1','var_data1'],
             ['rot_data2','rot_pc2','loD_data2','var_data2'],
             ['rot_data3','rot_pc3', 'loD_data3','var_data3']]
rot_styles= [['r.', 'r-'],
             ['g.', 'g-'],
             ['y.', 'y-']]

#Setup all layers
plotter.addLayer('orig_data')
plotter.addLayer('meta_data', kind='text')

for rot in rot_layers:
    for lay in rot:
        plotter.addLayer(lay)

plotter.arrangeFig([1,3], {'11':
                           {'name': 'BEFORE',
                            'scale': DOI,
                            'layers': ['orig_data','meta_data',
                                       'rot_data1', 'rot_pc1',
                                       'rot_data2', 'rot_pc2',
                                       'rot_data3', 'rot_pc3'],  # all layers will be selected
                            'axes_vars': ['x', 'y', 'z'],
                            'projection':'3d'},
                           '12':
                           {'name': 'AFTER',
                            'scale': [(-20,20),(-20,20)],
                            'layers': ['loD_data1',
                                       'loD_data2',
                                       'loD_data3'],  # all layers will be selected
                            'axes_vars': ['a', 'b']},
                           '13':
                           {'name': 'Variance by Components',
                            'scale': [(0,dim),(0,1)],
                            'layers': ['var_data1',
                                       'var_data2',
                                       'var_data3'],  # all layers will be selected
                            'axes_vars': ['x', 'y']},
                           })

gui.buildPlotter2D((8,8), with_times=False)

class ControlSys:
    def __init__(self, fig, data, rot_layers, rot_styles, d):
        self.fig = fig
        self.data = data
        self.rot_layers = rot_layers
        self.rot_styles = rot_styles
        self.d = d
        self.c = 0
        self.m = False
        self.fig.canvas.mpl_connect('key_press_event', self.keypress)

        print("Press left or right arrow keys to view different rotations of Hi-D data and their PC's.")
        print("Press m to display or hide all layers.")
        print("Press h to show or hide original data.")

    def keypress(self, event):
        if event.key == 'right':
            self.c += 1
        if event.key == 'left':
            self.c -= 1

        if event.key == 'left' or event.key == 'right':
            for rot in self.rot_layers:
                for lay in rot:
                    plotter.setLayer(lay, display= False)

            for i in range(len(self.rot_layers[0])):
                plotter.toggleDisplay(layer=self.rot_layers[self.c%len(self.rot_layers)][i]) #figure='Master',

        if event.key == 'm':
            self.m = not self.m
            for rot in self.rot_layers:
                for lay in rot:
                    plotter.setLayer(lay, figure='Master', display=self.m)

        if event.key == 'h':
            plotter.toggleDisplay(layer='orig_data', figure='Master')

        if event.key == 'up':
            self.d += 1
            loopPCA(self.data, self.d, self.rot_layers, self.rot_styles)

        if event.key == 'down':
            if self.d is not 2:
                self.d -= 1
                loopPCA(self.data, self.d, self.rot_layers, self.rot_styles)

        plotter.show(rebuild=False)

def disc():
    pts = sd.generate_ball(100, 2, 10)
    pts = np.concatenate((pts,np.zeros((100, 1))), axis=1)

    plotter.addData([pts[:,0], pts[:,1], pts[:,2]], layer='orig_data', style='b.')
    d=2

    ctrl_sys = ControlSys(gui.masterWin, pts, rot_layers, rot_styles, 2)

    loopPCA(pts, d, rot_layers, rot_styles)

def hypersphere(dim, rot_layers, rot_styles):
    pts = sd.generate_ball(100, dim, 10)
    plotter.addData([pts[:,0], pts[:,1], pts[:,2]], layer='orig_data', style='b.')
    d = 2

    ctrl_sys = ControlSys(gui.masterWin, pts, rot_layers, rot_styles, 2)

    loopPCA(pts, d, rot_layers, rot_styles)

hypersphere(6, rot_layers, rot_styles)
#disc()

halt= True