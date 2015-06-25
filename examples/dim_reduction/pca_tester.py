import pca_disc
from pca_disc import plotter, gui, compute, stretch, ControlSys, rotate_x, rotate_y, rotate_z, translate, noise
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

    X1 = noise(X1, 2, 0.3, 0, 10)
    X2 = noise(X2, 2, 0.3, 0, 10)
    X3 = noise(X3, 2, 0.3, 0, 10)

    X = [X1, X2, X3]

    rot_layers = ['rot1', 'rot2', 'rot3']
    rot_styles = ['r', 'g', 'y']

    pca_disc.setupDisplay(rot_layers, rot_styles, DOI)

    plotter.addData([pts[:,0], pts[:,1], pts[:,2]], layer='orig_data', style='b.')

    d=2
    ctrl_sys = ControlSys(gui.masterWin, X, rot_layers, rot_styles, 2)

    #for i in range(len(rot_layers)):
        #compute(X[i], d, rot_layers[i], rot_styles[i])

    return ctrl_sys

def hypersphere(dim):
    pts = sd.generate_ball(100, dim, 10)

    #Create and stretch different hypersphere "clusters":
    X1 = stretch(sd.generate_ball(133, dim, 10), 0, 1.2)
    X2 = sd.generate_ball(110, dim, 10)
    X3 = sd.generate_ball(95, dim, 10)

    X = [X1, X2, X3]

    clus_layers = ['clus1', 'clus2', 'clus3']
    clus_styles = ['r', 'g', 'y']

    pca_disc.setupDisplay(clus_layers, clus_styles, DOI)

    plotter.addData([pts[:,0], pts[:,1], pts[:,2]], layer='orig_data', style='b.')

    proj_vecsHI = pca_disc.ortho_proj_mat(len(X[0][0]), 3)
    proj_vecsLO = pca_disc.ortho_proj_mat(len(X[0][0]), 2)

    d = 2
    ctrl_sys = ControlSys(gui.masterWin, X, clus_layers, clus_styles, 2, proj_vecsLO, proj_vecsHI)

    return ctrl_sys

def iris():
    """
    Example for loading in real data from txt. Not yet working.
    """
    data = np.loadtxt('iris.data.txt', dtype= [('f0',float),('f1',float), ('f2',float),('f3',float), ('f4', object)], delimiter=',', unpack=False)
    X = {}

    for row in data:
        try:
            np.row_stack((X[row[4]], (np.array([row[i] for i in range(4)]))))
        except:
            X[row[4]] = (np.array([row[i] for i in range(4)]))

#ctrl_sys = disc()
ctrl_sys = hypersphere(6)

halt= True