from __future__ import absolute_import, print_function

from PyDSTool.Toolbox import data_analysis as da
from PyDSTool.Toolbox import synthetic_data as sd

import numpy as np

from fovea import *
from fovea.graphics import gui
from fovea.diagnostics import diagnostic_manager

import random

def translate(X, axis, amount):
    A = np.zeros((len(X), len(X[0])))
    A[:,axis] = np.ones(len(X))*amount
    Y = X + A
    return Y

def rotate_x(X, theta):
    R = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
    Y = np.dot(R, X.transpose())
    return Y.transpose()

def rotate_y(X, theta):
    R = np.array([[np.cos(theta), 0, -np.sin(theta)], [0, 1, 0], [np.sin(theta), 0, np.cos(theta)]])
    Y = np.dot(R, X.transpose())
    return Y.transpose()

def rotate_z(X, theta):
    R = np.array([[np.cos(theta), np.sin(theta), 0], [-np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    Y = np.dot(R, X.transpose())
    return Y.transpose()

def stretch(X, axis, amount):
    X[:,axis] = X[:,axis]*amount
    return X

def ortho_proj_mat(n, m):
    Z = npy.random.rand(n, m)
    Q, R = npy.linalg.qr(Z)
    return Q


plotter= gui.plotter

DOI = [(-10,10),(-10,10)]

plotter.clean() # in case rerun in same session
plotter.addFig('Master',
               title='PCA Disc',
               xlabel='x', ylabel='y',
               domain=DOI)

plotter.addLayer('orig_data')

plotter.addLayer('rot_data1')
plotter.addLayer('rot_pc1')

plotter.addLayer('rot_data2')
plotter.addLayer('rot_pc2')

plotter.addLayer('rot_data3')
plotter.addLayer('rot_pc3')

plotter.addLayer('loD_data')
plotter.addLayer('var_data')
plotter.addLayer('meta_data', kind='text')

plotter.arrangeFig([1,3], {'11':
                           {'name': 'Hi D',
                            'scale': DOI,
                            'layers': ['orig_data','meta_data',
                                       'rot_data1', 'rot_pc1',
                                       'rot_data2', 'rot_pc2',
                                       'rot_data3', 'rot_pc3'],  # all layers will be selected
                            'axes_vars': ['x', 'y', 'z'],
                            'projection':'3d'},
                           '12':
                           {'name': 'Lo D',
                            'scale': [(-10,10),(-10,10)],
                            'layers': 'loD_data',  # all layers will be selected
                            'axes_vars': ['x', 'y']},
                           '13':
                               {'name': 'Variance by Components',
                                'scale': [(0,10),(0,1)],
                                'layers': 'var_data',  # all layers will be selected
                                'axes_vars': ['x', 'y']},
                           })

gui.buildPlotter2D((8,8), with_times=False)

#Amount to translate data and along which axis.
trans_am = 15
trans_ax = 1

#Original data points. Set dim to an integer to view hyperspheres. Makes variance plot more interesting.
dim = 'disc'
#dim = 5

if dim == 'disc':
    pts = sd.generate_ball(100, 2, 10)
    pts = np.concatenate((pts,np.zeros((100, 1))), axis=1)
else:
    pts = sd.generate_ball(100, dim, 10)

pts = stretch(pts, 0, 1.7) #Stretching data causes one PC to capture a greater amount of variance.

plotter.addData([pts[:,0], pts[:,1], pts[:,2]], layer='orig_data', style='b.')

#Plot out several different rotations of the original data.
rot_layers= [['rot_data1','rot_pc1','r.', 'r-'], ['rot_data2','rot_pc2','g.', 'g-'], ['rot_data3','rot_pc3','y.', 'y-']]
for i in rot_layers:

    #If data high dimensional, create an arbitrary projection matrix so we can visualize.
    if(len(pts[0]) > 3):
        Y = pts
        Q3 = ortho_proj_mat(len(pts[0]), 3)

    else:
        Y = rotate_z(rotate_y(rotate_x(translate(pts, trans_ax, trans_am),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi))

    p = da.doPCA(Y, len(Y[0]), len(Y[0])) #Creates a pcaNode object.

    pcMat = p.get_projmatrix()

    pcPts = np.row_stack((pcMat.transpose()[0,:],np.zeros((1, len(pcMat))), pcMat.transpose()[1,:]))*10 #Format proj vecs as points.
    loPts = p._execute(Y) #Get dimensionality reduced data.

    if len(Y[0]) > 3:
        Y = npy.dot(Y, Q3)
        pcPts = npy.dot(pcPts, Q3)

    if len(loPts[0]) > 2:
        Q2 = ortho_proj_mat(len(loPts[0]), 2)
        loPts = npy.dot(loPts, Q2)

    #Create line plot for variance explained by each component.
    plotter.addData([range(len(p.d)), p.d/sum(p.d)], layer='var_data', style=i[3])

    #Create plot of high-dimensional data and its PC's.
    plotter.addData([Y[:,0], Y[:,1], Y[:,2]], layer=i[0], style=i[2])
    plotter.addData([pcPts[:,0], pcPts[:,1], pcPts[:,2]], layer=i[1], style=i[3])

    #Create plot of low-dimensional data.
    plotter.addData([loPts[:,0], loPts[:,1]], layer='loD_data', style=i[2])

    plotter.setLayer(i[0], figure='Master', display=False)
    plotter.setLayer(i[1], figure='Master', display=False)

plotter.auto_scale_domain(figure= 'Master')

plotter.show(rebuild=False)

c = 0
def keypress(event):
    global c
    c = c + 1

    if event.key == 'm':
        for i in rot_layers:
            plotter.setLayer(i[0], figure='Master', display= False)
            plotter.setLayer(i[1], figure='Master', display= False)

        plotter.toggleDisplay(layer=rot_layers[c%3][0], figure='Master')
        plotter.toggleDisplay(layer=rot_layers[c%3][1], figure='Master')

        plotter.show(rebuild=False)


gui.masterWin.canvas.mpl_connect('key_press_event', keypress)

print("Press 'm' to view different rotations of Hi-D data and their PC's")
halt=True