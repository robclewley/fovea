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

def noise(X, axis, percent, loc, scale):
    for i in range(0, round(len(X)*percent)):
        X[i][axis] = np.random.normal(loc, scale, 1)
    return X

def ortho_proj_mat(n, m):
    Z = npy.random.rand(n, m)
    Q, R = npy.linalg.qr(Z)
    return Q

print("Press left or right arrow keys to view different rotations of Hi-D data and their PC's.")
print("Press m to display or hide all layers.")
print("Press h to show or hide original data.")

plotter= gui.plotter

DOI = [(-10,10),(-10,10)]

plotter.clean() # in case rerun in same session
plotter.addFig('Master',
               title='PCA Disc',
               xlabel='x', ylabel='y',
               domain=DOI)

#Original data points. Set dim to an integer to view hyperspheres. Makes variance plot more interesting.
#dim = 'disc'
dim = 6

if dim == 2:
    pts = sd.generate_ball(100, 2, 10)
    pts = np.concatenate((pts,np.zeros((100, 1))), axis=1)
else:
    pts = sd.generate_ball(100, dim, 10)

#Setup all layers
plotter.addLayer('orig_data')

plotter.addLayer('rot_data1')
plotter.addLayer('rot_pc1')
plotter.addLayer('loD_data1')
plotter.addLayer('var_data1')

plotter.addLayer('rot_data2')
plotter.addLayer('rot_pc2')
plotter.addLayer('loD_data2')
plotter.addLayer('var_data2')

plotter.addLayer('rot_data3')
plotter.addLayer('rot_pc3')
plotter.addLayer('loD_data3')
plotter.addLayer('var_data3')

plotter.addLayer('meta_data', kind='text')

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
                            'axes_vars': ['x', 'y']},
                           '13':
                               {'name': 'Variance by Components',
                                'scale': [(0,dim),(0,1)],
                                'layers': ['var_data1',
                                           'var_data2',
                                           'var_data3'],  # all layers will be selected
                                'axes_vars': ['x', 'y']},
                           })

gui.buildPlotter2D((8,8), with_times=False)

#Amount to translate data and along which axis.
trans_am = 15
trans_ax = 1

pts = stretch(pts, 0, 1.2) #Stretching data causes one PC to capture a greater amount of variance.
pts = stretch(pts, 1, 1.2)
pts = stretch(pts, 2, 1.2)
#pts = noise(pts, 2, 0.8, 0, 1.5)

plotter.addData([pts[:,0], pts[:,1], pts[:,2]], layer='orig_data', style='b.')

#Plot out several different rotations of the original data.
rot_layers= [['rot_data1','rot_pc1','loD_data1','var_data1'],
             ['rot_data2','rot_pc2','loD_data2','var_data2'],
             ['rot_data3','rot_pc3', 'loD_data3','var_data3']]
rot_styles= [['r.', 'r-'],
             ['g.', 'g-'],
             ['y.', 'y-']]

def loopPCA(pts, new_dim, rot_layers, rot_styles):
    for i in range(len(rot_layers)):

        for j in rot_layers[i]:
            plotter.setLayer(j, figure='Master', data={})

        #If data high dimensional, create an arbitrary projection matrix so we can visualize.
        if(len(pts[0]) > 3):
            Y = pts
            Q3 = ortho_proj_mat(len(pts[0]), 3)

        else:
            Y = rotate_z(rotate_y(rotate_x(translate(pts, trans_ax, trans_am),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi)),random.uniform(0, 2*np.pi))

        p = da.doPCA(Y, len(Y[0]), len(Y[0])) #Creates a pcaNode object.

        pcMat = p.get_projmatrix()

        pcPts = np.concatenate(([pcMat.transpose()[0,:]*15],
                        [-pcMat.transpose()[0,:]*15]), axis=0)
        for j in range(1, new_dim):
            pcPts = np.concatenate((pcPts, [pcMat.transpose()[j,:]*15],
                            [-pcMat.transpose()[j,:]*15]), axis=0)

        loPts = p._execute(Y, new_dim) #Get dimensionality reduced data.

        if len(Y[0]) > 3:
            Y = npy.dot(Y, Q3)
            pcPts = npy.dot(pcPts, Q3) #Will this work with new pcPts?

        if len(loPts[0]) > 2:
            Q2 = ortho_proj_mat(len(loPts[0]), 2)
            loPts = npy.dot(loPts, Q2)

        #Create line plot for variance explained by each component.
        plotter.addData([range(len(p.d)), p.d/sum(p.d)], layer=rot_layers[i][3], style=rot_styles[i][1]+"o")

        #Create plot of high-dimensional data and its PC's.
        plotter.addData([Y[:,0], Y[:,1], Y[:,2]], layer= rot_layers[i][0], style=rot_styles[i][0])

        for j in range(0, len(pcPts), 2):
            plotter.addData([pcPts[j+0:j+2,0], pcPts[j+0:j+2,1], pcPts[j+0:j+2,2]], layer= rot_layers[i][1], style= rot_styles[i][1])
            plotter.addData([pcPts[j+2:j+4,0], pcPts[j+2:j+4,1], pcPts[j+2:j+4,2]], layer= rot_layers[i][1], style= rot_styles[i][1])

        #Create plot of low-dimensional data.
        plotter.addData([loPts[:,0], loPts[:,1]], layer=rot_layers[i][2], style=rot_styles[i][0])

        for j in rot_layers[i]:
            plotter.setLayer(j, figure='Master', display=False)

    print("Variance Explained:")
    print(sum(p.d[0:new_dim])/sum(p.d))

    plotter.show(rebuild=False)

c = 0
m = False
d = 2

loopPCA(pts, 2, rot_layers, rot_styles)

def keypress(event):
    global c
    global m
    global d

    if event.key == 'right':
        c = c + 1
    if event.key == 'left':
        c = c - 1

    if event.key == 'left' or event.key == 'right':
        for rot in rot_layers:
            for lay in rot:
                plotter.setLayer(lay, figure='Master', display= False)

        for i in range(len(rot_layers[0])):
            plotter.toggleDisplay(layer=rot_layers[c%len(rot_layers)][i], figure='Master')

    if event.key == 'm':
        m = not m
        for rot in rot_layers:
            for lay in rot:
                plotter.setLayer(lay, figure='Master', display=m)

    if event.key == 'h':
        plotter.toggleDisplay(layer='orig_data', figure='Master')

    if event.key == 'up':
        d = d + 1
        plotter.addVLine( d, layer='var_data1')
        loopPCA(pts, d, rot_layers, rot_styles)

    if event.key == 'down':
        if d is not 2:
            d = d - 1
            loopPCA(pts, d, rot_layers, rot_styles)

    plotter.show(rebuild=False)

gui.masterWin.canvas.mpl_connect('key_press_event', keypress)

halt=True
