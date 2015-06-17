from pca_disc import plotter, gui, loopPCA, stretch, setupPCAlayers, ControlSys
from fovea import *

from PyDSTool.Toolbox import synthetic_data as sd

import numpy as np
import __future__

#Plot out several different rotations of the original data.
rot_layers= [['rot_data1','rot_pc1','loD_data1','var_data1'],
             ['rot_data2','rot_pc2','loD_data2','var_data2'],
             ['rot_data3','rot_pc3', 'loD_data3','var_data3']]
rot_styles= [['r.', 'r-'],
             ['g.', 'g-'],
             ['y.', 'y-']]

DOI = [(-10,10),(-10,10)]

setupPCAlayers(rot_layers, rot_styles, DOI)

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