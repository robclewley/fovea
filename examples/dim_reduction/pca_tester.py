import pca_disc
from pca_disc import plotter, gui, loopPCA, stretch, ControlSys, rotate_x, rotate_y, rotate_z, translate
from fovea import *

from PyDSTool.Toolbox import synthetic_data as sd

import random
import numpy as np
import __future__

DOI = [(-10,10),(-10,10)]

def disc():
    pts = sd.generate_ball(100, 2, 10)
    pts = np.concatenate((pts,np.zeros((100, 1))), axis=1)

    trans_am = 15
    trans_ax = 1

    #Produce rotated clusters:
    X1 = rotate_z(rotate_y(rotate_x(translate(pts, trans_ax, trans_am),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi))
    X2 = rotate_z(rotate_y(rotate_x(translate(pts, trans_ax, trans_am),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi))
    X3 = rotate_z(rotate_y(rotate_x(translate(pts, trans_ax, trans_am),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi))

    X = [X1, X2, X3]

    #Plot out several different rotations of the original data.
    rot_layers= [['rot_data1','rot_pc1','loD_data1','var_data1'],
                 ['rot_data2','rot_pc2','loD_data2','var_data2'],
                 ['rot_data3','rot_pc3', 'loD_data3','var_data3']]
    rot_styles= [['r.', 'r-'],
                 ['g.', 'g-'],
                 ['y.', 'y-']]

    pca_disc.setupPCAlayers(rot_layers, rot_styles, DOI)

    plotter.addData([pts[:,0], pts[:,1], pts[:,2]], layer='orig_data', style='b.')

    d=2
    ctrl_sys = ControlSys(gui.masterWin, pts, rot_layers, rot_styles, 2)

    for i in range(len(rot_layers)):
        loopPCA(X[i], d, rot_layers[i], rot_styles[i])

def hypersphere(dim, clus_layers, clus_styles):
    pts = sd.generate_ball(100, dim, 10)

    #Produce rotated clusters:
    X1 = stretch(sd.generate_ball(133, dim, 10), 0, 1.2)
    X2 = sd.generate_ball(110, dim, 10)
    X3 = sd.generate_ball(95, dim, 10)

    X = [X1, X2, X3]

    #Plot out several different hyperspheres.
    clus_layers= [['clus1_data','clus1_pc','clus1_loD','clus1_var'],
                 ['clus2_data','clus2_pc','clus2_loD','clus2_var'],
                 ['clus3_data','clus3_pc','clus3_loD','clus3_var']]
    clus_styles= [['r.', 'r-'],
                 ['g.', 'g-'],
                 ['y.', 'y-']]


    pca_disc.setupPCAlayers(clus_layers, clus_styles, DOI)

    plotter.addData([pts[:,0], pts[:,1], pts[:,2]], layer='orig_data', style='b.')

    d = 4
    ctrl_sys = ControlSys(gui.masterWin, X, clus_layers, clus_styles, 2)

    proj_vecsHI = pca_disc.ortho_proj_mat(len(X[0][0]), 3)
    proj_vecsLO = pca_disc.ortho_proj_mat(d, 2)

    for i in range(len(clus_layers)):
        loopPCA(X[i], d, clus_layers[i], clus_styles[i], proj_vecsLO, proj_vecsHI)

#hypersphere(6, [], [])
disc()

halt= True