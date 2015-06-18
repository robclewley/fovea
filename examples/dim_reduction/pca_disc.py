from __future__ import absolute_import, print_function

from PyDSTool.Toolbox import data_analysis as da

import numpy as np

from fovea import *
from fovea.graphics import gui
from fovea.diagnostics import diagnostic_manager

import random

plotter= gui.plotter

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

def loopPCA(X, new_dim, layers, styles, proj_vecsLO = None, proj_vecsHI = None):

    p = da.doPCA(X, len(X[0]), len(X[0])) #Creates a pcaNode object.

    pcMat = p.get_projmatrix()

    pcPts = np.concatenate(([pcMat.transpose()[0,:]*15],
                            [-pcMat.transpose()[0,:]*15]), axis=0)

    for j in range(1, new_dim):
        pcPts = np.concatenate((pcPts, [pcMat.transpose()[j,:]*15],
                                [-pcMat.transpose()[j,:]*15]), axis=0)

    Y = p._execute(X, new_dim) #Get dimensionality reduced data.

    #Produce proj vecs if none supplied
    if proj_vecsHI is None:
        proj_vecsHI = ortho_proj_mat(len(X[0]), 3)


    if proj_vecsLO is None:
        proj_vecsLO = ortho_proj_mat(len(Y[0]), 2)

    for layer in layers:
        plotter.setLayer(layer, figure='Master', data={})

    #If data are high-dimensional, use the projection vectors.
    if len(X[0]) > 3:
        X = npy.dot(X, proj_vecsHI)
        pcPts = npy.dot(pcPts, proj_vecsHI)

    if len(Y[0]) > 2:
        proj_vecsLO = proj_vecsLO[0:new_dim, :] #Scrape the needed rows from projection matrix, so we can page between dimensionalities.
        Y = npy.dot(Y, proj_vecsLO)

    #Create line plot for variance explained by each component.
    plotter.addData([range(len(p.d)), p.d/sum(p.d)], layer=layers[3], style=styles[1]+"o")

    #Create plot of high-dimensional data and its PC's.
    plotter.addData([X[:,0], X[:,1], X[:,2]], layer= layers[0], style=styles[0])

    for j in range(0, len(pcPts), 2):
        plotter.addData([pcPts[j+0:j+2,0], pcPts[j+0:j+2,1], pcPts[j+0:j+2,2]], layer= layers[1], style= styles[1])
        plotter.addData([pcPts[j+2:j+4,0], pcPts[j+2:j+4,1], pcPts[j+2:j+4,2]], layer= layers[1], style= styles[1])

    #Create plot of low-dimensional data.
    plotter.addData([Y[:,0], Y[:,1]], layer=layers[2], style=styles[0])

    for layer in layers:
        plotter.setLayer(layer, figure='Master', display=False)

    print("Variance Explained:")
    print(sum(p.d[0:new_dim])/sum(p.d))

    plotter.show(rebuild=False)


def setupPCAlayers(rot_layers, rot_styles, DOI):
    plotter.clean() # in case rerun in same session
    plotter.addFig('Master',
                   title='PCA Disc',
                   xlabel='x', ylabel='y',
                   domain=DOI)

    #Setup all layers
    plotter.addLayer('orig_data')
    plotter.addLayer('meta_data', kind='text')

    for rot in rot_layers:
        for lay in rot:
            plotter.addLayer(lay)

    plotter.arrangeFig([1,3], {'11':
                               {'name': 'BEFORE',
                                'scale': DOI,
                                'layers': sum([rot_layers[i][0:2] for i in range(3)],['orig_data']),  # all layers will be selected
                                'axes_vars': ['x', 'y', 'z'],
                                'projection':'3d'},
                               '12':
                               {'name': 'AFTER',
                                'scale': [(-20,20),(-20,20)],
                                'layers': [rot_layers[i][2] for i in range(3)],  # all layers will be selected
                                'axes_vars': ['a', 'b']},
                               '13':
                               {'name': 'Variance by Components',
                                'scale': [(0,10),(0,1)],
                                'layers': [rot_layers[i][3] for i in range(3)],  # all layers will be selected
                                'axes_vars': ['x', 'y']},
                               })

    gui.buildPlotter2D((8,8), with_times=False)

class ControlSys:
    def __init__(self, fig, data, clus_layers, clus_styles, d, proj_vecsLO=None, proj_vecsHI=None):
        self.fig = fig
        self.data = data
        self.clus_layers = clus_layers
        self.clus_styles = clus_styles
        self.proj_vecsLO = proj_vecsLO
        self.proj_vecsHI = proj_vecsHI
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

            for clus in self.clus_layers:
                for lay in clus:
                    plotter.setLayer(lay, display= False)

            for i in range(len(self.clus_layers[0])):
                plotter.toggleDisplay(layer=self.clus_layers[self.c%len(self.clus_layers)][i]) #figure='Master',

        if event.key == 'm':
            self.m = not self.m
            for clus in self.clus_layers:
                for lay in clus:
                    plotter.setLayer(lay, figure='Master', display=self.m)

        if event.key == 'h':
            plotter.toggleDisplay(layer='orig_data', figure='Master')

        if event.key == 'up':
            self.d += 1
            for i in range(len(self.clus_layers)):
                loopPCA(self.data[i], self.d, self.clus_layers[i], self.clus_styles[i], self.proj_vecsLO, self.proj_vecsHI) #Issue: Going to be randomly generating proj_vecs.

        if event.key == 'down':
            if self.d is not 2:
                self.d -= 1
                for i in range(len(self.clus_layers)):
                    loopPCA(self.data[i], self.d, self.clus_layers[i], self.clus_styles[i], self.proj_vecsLO, self.proj_vecsHI)

        plotter.show(rebuild=False)

halt=True
