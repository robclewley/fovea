"""
This is the main run script for simple demo involving
Hodgkin-Huxley analysis.
"""

from PyDSTool.Toolbox.dssrt import *
from PyDSTool.Toolbox.phaseplane import *
import PyDSTool as dst

import matplotlib.pyplot as plt
import sys

from fovea.graphics import gui
from model_config import man as modelManager
from model_config import name, header
from HH_neuron import getHH_DSSRT, computePPlaneObjects, do_traj
from fovea.common import castNull, castNullArray
from math import *

from scipy.optimize import fsolve

## ----- ----- ----- ----- ----- ----- ##
## BUILD GENERATOR OBJECT              ##
## ----- ----- ----- ----- ----- ----- ##

model = modelManager.instances[name]
gen = list(model.registry.values())[0]

# global define for convenience
plotter = gui.plotter

# ------------------------

def clip_to_pt():
    """Extract clipboard point from gui to a dictionary"""
    pt = dst.filteredDict(gui.capturedPts['Master'], ['V', 'm', 'n'])
    return {'V': pt['V'], 'Na.m': pt['m'], 'K.n': pt['n']}


class PPcallback(object):
    """
    Dynamic figure axes class to support state-dependent user-interactive callbacks
    """
    def __init__(self, xvar, num_x_points=30, num_y_points=30,
                 nullcX_style=None, nullcY_style=None,
                 vel_arrow_scale=1):
        self.nully = None
        self.nullx = None
        self.num_x_points = num_x_points
        self.num_y_points = num_y_points
        if nullcX_style is None:
            self.nullcX_style = 'b-'
        else:
            self.nullcX_style = nullcX_style
        if nullcY_style is None:
            self.nullcY_style = 'r-'
        else:
            self.nullcY_style = nullcY_style
        self.last_scale = None
        self.vel_arrow_scale = vel_arrow_scale
        self.first_call = True  # is reset by __call__


    def dQ_dt(self, Qstr, ix, points):
        """
        Utility to find finite difference of any named quantity in given points
        at index ix
        """
        if ix == 0:
            ix = 1
        pt1 = points[ix]
        pt0 = points[ix-1]
        t1 = points.indepvararray[ix]
        t0 = points.indepvararray[ix-1]
        return (pt1[Qstr]-pt0[Qstr])/(t1-t0)


class PPcallback_m(PPcallback):
    def __call__(self, time, hard_reset=False):
        """Callback 'function' to take care of refreshing and re-computing
        phase plane sub-plot when time is changed interactively.
        """
        #print("\n\npplane call back, mouseUp =", gui._mouseUp)

        fig_struct = plotter.figs['Master']
        # combine all layer information
        dynamicData = fig_struct.layers['nullclines_mV'].data.copy()
        dynamicData.update(fig_struct.layers['horiz_PP'].data)

        ax = gui.dynamicPlots['nullclines_mV']
        #sc = fig_struct.layers['nullclines'].scale
        sc = [ax.get_xlim(), ax.get_ylim()]

        if hard_reset:
            force = True
            preComputed = False
            print("\n HARD REFRESH for phase plane")
        else:
            preComputed = False
            force = False
            # REPLACE WITH A PROPER CACHE STRUCTURE
            for key, val in dynamicData.items():
                # dynamicData.keys are 'yNull_<time>' or 'xNull' or keys from horiz_PP etc.
                # re-use computed nullclines if time is in "cache", i.e. it shows up in the keys
                if key[6:] == str(time):
                    # cache hit!
                    val['display'] = True  # not clear if this updates original data structure after copy
                    # also check to see whether has been rescaled
                    if self.last_scale == sc:
                        preComputed = True
                    else:
                        force = True
                elif key[:5] != 'xNull':
                    # Use != to clean up collision lines and ohter nullcline that are not for
                    # this time value.
                    # switch off y-nullcline (V_inf) display for the other times
                    # yNull stays constant so keep that display=True
                    val['display'] = False

        pt = gui.points[gui.ix]

        p = fig_struct.layers['points_mV']
        p.display = True

        dV_dt = (pt['vinf']-pt['V'])/pt['tauv']
        dm_dt = (pt['Na.minf']-pt['Na.m'])/pt['Na.taum']
        dn_dt = (pt['K.ninf']-pt['K.n'])/pt['K.taun']
        gui.add_data_points([[pt['Na.m'], pt['Na.m']+dm_dt*self.vel_arrow_scale],
                         [pt['V'], pt['V']+dV_dt*self.vel_arrow_scale]],
                        layer='state_vel_mV', name='state', style=vel_vec_style, force=True)

        if self.first_call:
            gui.add_data_points([gui.points['Na.m'], gui.points['V']],
                            layer='vfp_mV', name='traj', style='y')

        # Virtual fixed point and linearized nullclines
        if 'fast_m' in model.name:
            with_jac = False
            do_fps = False
            fast_vars = ['Na.m']
        else:
            with_jac = False
            do_fps = False
            fast_vars = None

        # update (or create) points
        try:
            gui.set_point('state_pt', Point2D(pt['Na.m'], pt['V']), 'points_mV')
            gui.set_point('vinf_pt', Point2D(pt['Na.m'], pt['vinf']), 'points_mV')
        except KeyError:
            gui.add_data_points(Point2D(pt['Na.m'], pt['V']), coorddict = {'x':
                         {'y':'y', 'style':'ko', 'layer':'points_mV', 'name':'state_pt'}})
            gui.add_data_points(Point2D(pt['Na.m'], pt['vinf']),coorddict = {'x':
                         {'y':'y', 'style':'bx', 'layer':'points_mV', 'name':'vinf_pt'}})

        d = fig_struct.layers['nullclines_mV'].data

        if not preComputed and gui._mouseUp:
##            print("\nComputing phase plane...")
##            print("  Current time = %.4f" % (time))

            if self.nullx is None or force:
                # compute m nullcline this once
                only_var = None
            else:
                only_var = 'V'

            # refresh wait notification
            ax.text(0.05, 0.95, 'wait', transform=ax.transAxes, fontsize=22,
                    color='r', fontweight='bold', va='top')
            gui.masterWin.canvas.draw()

            # comment out for testing - use surrogate below
            nulls = computePPlaneObjects(gen, 'Na.m', 'V', state=pt,
                                         num_x_points=self.num_x_points,
                                         num_y_points=self.num_y_points,
                                         only_var=only_var, with_jac=with_jac,
                                         do_fps=do_fps, fast_vars=fast_vars,
                                         subdomain={'V': sc[1],
                                                    'Na.m': sc[0]})

            # Surrogate data - much faster to test with
            #self.nully = [[-100+time, -50+time/10., 0], [0.1, 0.4, 0.8]]
            #self.nullx = [[-130, -80, 50], [0.2, 0.3, 0.4]]

            self.nully = castNullArray(nulls['nullcY'])
            gui.add_data_points(self.nully, layer='nullclines_mV', style=self.nullcY_style,
                            name='yNull_'+str(time), force=force)

            # delete update 'wait' notice
            ax.texts = []
            #ax.clear()
            gui.clear_axes(ax)

            if only_var is None:
                # nullx is added second so will be the second line
                self.nullx = castNullArray(nulls['nullcX'])
                gui.add_data_points(self.nullx, layer='nullclines_mV',
                                style=self.nullcX_style,
                                name='xNull', force=force)

            #if force:
            #    rescale = sc
            #else:
            #    rescale = None
            gui.build_layers(['nullclines_mV', 'horiz_PP',  'points_mV',
                             'state_vel_mV', 'vfp_mV'],
                            ax, rescale=sc, figure='Master')

            self.last_scale = sc
##            print("  Phase plane rebuild completed.\n")
        else:
            # just refresh display with the current selected data
            gui.clear_axes(ax)
            #if force:
            #    rescale = sc
            #else:
            #    rescale = None
            gui.build_layers(['nullclines_mV', 'horiz_PP',  'points_mV',
                             'state_vel_mV', 'vfp_mV'],
                            ax, rescale=sc, figure='Master')
            self.last_scale = sc

        gui.masterWin.canvas.draw()
        self.first_call = False


class PPcallback_n(PPcallback):
    def __call__(self, time, hard_reset=False):
        """Callback 'function' to take care of refreshing and re-computing
        phase plane sub-plot when time is changed interactively.
        """
        #print("\n\npplane call back, mouseUp =", gui._mouseUp)

        fig_struct = plotter.figs['Master']
        # combine all layer information
        dynamicData = fig_struct.layers['nullclines_nV'].data.copy()

        ax = gui.dynamicPlots['nullclines_nV']
        #sc = fig_struct.layers['nullclines'].scale
        sc = [ax.get_xlim(), ax.get_ylim()]

        if hard_reset:
            force = True
            preComputed = False
            print("\n HARD REFRESH for phase plane")
        else:
            preComputed = False
            force = False
            # REPLACE WITH A PROPER CACHE STRUCTURE
            for key, val in dynamicData.items():
                # dynamicData.keys are 'yNull_<time>' or 'xNull' or keys from horiz_PP etc.
                # re-use computed nullclines if time is in "cache", i.e. it shows up in the keys
                if key[6:] == str(time):
                    # cache hit!
                    val['display'] = True  # not clear if this updates original data structure after copy
                    # also check to see whether has been rescaled
                    if self.last_scale == sc:
                        preComputed = True
                    else:
                        force = True
                elif key[:5] != 'xNull':
                    # Use != to clean up collision lines and ohter nullcline that are not for
                    # this time value.
                    # switch off y-nullcline (V_inf) display for the other times
                    # yNull stays constant so keep that display=True
                    val['display'] = False

        pt = gui.points[gui.ix]

        p = fig_struct.layers['points_nV']
        p.display = True

        dV_dt = (pt['vinf']-pt['V'])/pt['tauv']
        dm_dt = (pt['Na.minf']-pt['Na.m'])/pt['Na.taum']
        dn_dt = (pt['K.ninf']-pt['K.n'])/pt['K.taun']
        gui.add_data_points([[pt['K.n'], pt['K.n']+dn_dt*self.vel_arrow_scale],
                         [pt['V'], pt['V']+dV_dt*self.vel_arrow_scale]],
                        layer='state_vel_nV', name='state', style=vel_vec_style, force=True)

        if self.first_call:
            gui.add_data_points([gui.points['K.n'], gui.points['V']],
                            layer='vfp_nV', name='traj', style='y')
            gui.add_data_points([gui.points['K.n'], gui.points['vinf']],
                            layer='vfp_nV', name='quasiVnull', style='m--')
##            vs = np.linspace(sc[1][0], sc[1][1], 50)
##            x = dict(pt).copy()
##
##            def vinf(n, v):
##                x['K.n'] = n
##                x['V'] = v
##                x['Na.m'] = gen.auxfns.Na_dssrt_fn_minf(v)
##                # assume autonomous system
##                return model.Rhs(0, x, asarray=False)['V']
##
##            vinfs_inv_n = [fsolve(vinf, gen.auxfns.K_dssrt_fn_ninf(v), args=(v,)) for v in vs]
##            plotter.add_data([vinfs_inv_n, vs], layer='vfp_nV', name='vinf_fastm', style='b--')

        # Virtual fixed point and linearized nullclines
        if 'fast_m' in model.name:
            with_jac = False
            do_fps = False
            fast_vars = ['Na.m']
        else:
            with_jac = False
            do_fps = False
            fast_vars = None

        # update (or create) points
        try:
            gui.set_point('state_pt', Point2D(pt['K.n'], pt['V']), 'points_nV')
            gui.set_point('vinf_pt', Point2D(pt['K.n'], pt['vinf']), 'points_nV')
        except KeyError:
            gui.add_data_points(Point2D(pt['K.n'], pt['V']),coorddict = {'x':
                {'y':'y', 'style':'ko', 'layer':'points_nV', 'name':'state_pt'}})
            gui.add_data_points(Point2D(pt['K.n'], pt['vinf']), coorddict = {'x':
                {'y':'y', 'style':'bx', 'layer':'points_nV', 'name':'vinf_pt'}})

        d = fig_struct.layers['nullclines_nV'].data

        if not preComputed and gui._mouseUp:
##            print("\nComputing phase plane...")
##            print("  Current time = %.4f" % (time))

            if self.nullx is None or force:
                # compute m nullcline this once
                only_var = None
            else:
                only_var = 'V'

            # refresh wait notification
            ax.text(0.05, 0.95, 'wait', transform=ax.transAxes, fontsize=22,
                    color='r', fontweight='bold', va='top')
            gui.masterWin.canvas.draw()

            # comment out for testing - use surrogate below
            nulls = computePPlaneObjects(gen, 'K.n', 'V', state=pt,
                                         num_x_points=self.num_x_points,
                                         num_y_points=self.num_y_points,
                                         only_var=only_var, with_jac=with_jac,
                                         do_fps=do_fps, fast_vars=fast_vars,
                                         subdomain={'V': sc[1],
                                                    'K.n': sc[0]})

            # Surrogate data - much faster to test with
            #self.nully = [[-100+time, -50+time/10., 0], [0.1, 0.4, 0.8]]
            #self.nullx = [[-130, -80, 50], [0.2, 0.3, 0.4]]

            self.nully = castNullArray(nulls['nullcY'])
            gui.add_data_points(self.nully, layer='nullclines_nV', style=self.nullcY_style,
                              name='yNull_'+str(time), force=force)

            # delete update 'wait' notice
            ax.texts = []
            #ax.clear()
            gui.clear_axes(ax)

            if only_var is None:
                # nullx is added second so will be the second line
                self.nullx = castNullArray(nulls['nullcX'])
                gui.add_data_points(self.nullx, layer='nullclines_nV',
                                  style=self.nullcX_style, name='xNull', force=force)

            #if force:
            #    rescale = sc
            #else:
            #    rescale = None
            gui.build_layers(['nullclines_nV', 'points_nV', 'state_vel_nV', 'vfp_nV'],
                            ax, rescale=sc, figure='Master')

            self.last_scale = sc
##            print("  Phase plane rebuild completed.\n")
        else:
            # just refresh display with the current selected data
            gui.clear_axes(ax)
            #if force:
            #    rescale = sc
            #else:
            #    rescale = None
            gui.build_layers(['nullclines_nV', 'points_nV', 'state_vel_nV', 'vfp_nV'],
                            ax, rescale=sc, figure='Master')
            self.last_scale = sc

        gui.masterWin.canvas.draw()
        self.first_call = False


# ------------------------------------------------------------------------

## Set dssrt_name to be for saved DSSRT data file
## Change for different parameter sets or just default to model name

dssrt_name = name

if 'typeI' in name and 'typeII' not in name:
    ### FOR TYPE I H-H ONLY
    Kgmax = 100 # 100 # 80 original
    dssrt_name = name+'_gK%i' % Kgmax
    model.set(pars={'K.g': Kgmax,
                    'Na.g': 50})
else:
    ### FOR TYPE II H-H ONLY
    Kgmax = 36 #39 or 42 with fast m # 36 original
    dssrt_name = name+'_gK%i' % Kgmax
    model.set(pars={'K.g': Kgmax,
                    'Ib.Ibias': 8.})

dV = 0.2

##test_ic = {'K.n': 0.37220277852490802,
##           'Na.m': 0.080387043479386036,
##           'V': -59.5}
##model.set(ics=test_ic)


## ----- ----- ----- ----- ----- ----- ##
## GET GENERATOR TRAJECTORY            ##
## ----- ----- ----- ----- ----- ----- ##

orig_ics = model.query('ics')

if 'no_h' in name:
    # no periodic orbit, just simulate for 12 ms
    if 'typeI' in name and 'typeII' not in name:
        t_end = 9
    else:
        t_end = 12
    model.set(tdata=[0, t_end])
    model.compute('ref')
    ref_traj = model['ref']
else:
    # get periodic orbit
    t_end = 20
    ref_traj, ref_pts, ref_tmin, ref_tmax = do_traj(model, t_end,
                                                    do_plot=False)

# re-sample traj at constant dt and declare to GUI
#trajPts = ref_traj.sample(dt=0.01)[:-40] #  cheap way to avoid overlap from pts not being periodic
trajPts = ref_traj.sample(dt=0.01)[:len(ref_traj.sample(dt=0.01))-40] #[:-40] syntax not working in python 3
gui.add_time_from_points(trajPts)


## ----- ----- ----- ----- ----- ----- ##
## CREATE DIAGNOSTIC OBJECT            ##
## ----- ----- ----- ----- ----- ----- ##

gui.clean()
gui.add_fig('Master', title='Geometric Dynamic Analysis: '+dssrt_name,
               tdom=[0, t_end], domain=[(-100,50), (0,1)])
coorddict = {'V':
             {'x':'t', 'layer':'V','name':'V', 'style':'k-'},
             'vinf':
             {'x':'t', 'layer':'V','name':'Vinf', 'style':'k:'},
             'Na.m':
             {'x':'t', 'layer':'activs', 'name':'m', 'style':'g--'},
             'Na.minf':
             {'x':'t', 'layer':'activs', 'name':'minf', 'style':'g--'},
             'K.n':
             {'x':'t', 'layer':'activs', 'name':'n', 'style':'r-'},
             'K.ninf':
             {'x':'t', 'layer':'activs', 'name':'ninf', 'style':'r--'},
             'tauv':
             {'x':'t','layer':'activs','name':'tauv', 'style':'b:'},
             'Na.taum':
             {'x':'t', 'layer':'activs','name':'taum', 'style':'g:'},
             'K.taun':
             {'x':'t', 'layer':'activs','name':'taun', 'style':'r:'}
             }
gui.add_data_points(trajPts, coorddict = coorddict)

print("Key for activations / time scales window")
print("  Activations: line=activation, dashed=asymptotic")
print("  Time scales: dots")
print("Na: green")
print("K: red")
print("V: blue")

## ----- ----- ----- ----- ----- ----- ##
## COMPUTE V-m PHASE PLANE             ##
## ----- ----- ----- ----- ----- ----- ##

# start at t = 2ms
gui.set_time(2)

# global style defs
vel_vec_style = {'color': 'k', 'linewidth': 2, 'linestyle': '-'}
vinf_vec_style = {'color': 'b', 'linewidth': 2, 'linestyle': '-'}
horiz_style = {'color': 'k', 'linestyle': '--', 'linewidth': 2}

def make_layer(xvar):
    if xvar == 'Na.m':
        suffix = 'mV'
    else:
        suffix = 'nV'
    PP_layer_name = 'nullclines_'+suffix
    gui.add_layer(PP_layer_name, dynamic=True)
    if xvar == 'Na.m':
        gui.add_layer('horiz_PP')
        nullcX_style = 'g-'
        PPclass = PPcallback_m
    else:
        # no horizon layer for K.n
        nullcX_style = 'r-'
        PPclass = PPcallback_n
    PPplot = PPclass(xvar, nullcY_style = {'color': 'b', 'linestyle': '-', 'linewidth': 1},
                                                      nullcX_style=nullcX_style)
    gui.dynamicPlotFns[PP_layer_name] = PPplot
    gui.add_layer('points_'+suffix)
    gui.add_layer('state_vel_'+suffix)
    gui.add_layer('vfp_'+suffix)

make_layer('Na.m')
make_layer('K.n')

# sub-plot specs: (row, col) integer coords start at top left
dPlot11 = {'name': 'Trajectory',
           'layers': ['V'],
           'scale': [None, [-85, 45]],
           'axes_vars': ['t', 'V'] }

dPlot12 = {'name': 'Activations, Time scales',
           'layers': ['activs'],
           'scale': [None, [0,1]],
           'axes_vars': ['t', 'no units, ms'] }

pp1_name = 'Na-V Phaseplane'
pp1_dom = [[0,1], [-75,50]]
pp1_vars = ['m', 'V']
pp2_name = 'K-V Phaseplane'
pp2_dom = [[0.,1], [-75,50]]
pp2_vars = ['n', 'V']

# ISSUE: Rename 'scale' to 'domain' or 'extent'
dPlot21 = {'name': pp1_name,
           'scale': pp1_dom,
           'layers': ['nullclines_mV', 'horiz_PP', 'points_mV', 'state_vel_mV', 'vfp_mV'],
           'axes_vars': pp1_vars}

dPlot22 = {'name': pp2_name,
           'scale': pp2_dom,
           'layers': ['nullclines_nV', 'points_nV', 'state_vel_nV', 'vfp_nV'],
           'axes_vars': pp2_vars}


dPlot_dict = {'11': dPlot11, '12': dPlot12, '21': dPlot21, '22': dPlot22}

gui.setup(dPlot_dict, size=(14, 8))

gui.show_legends(subplot='Times')

gui.show()

gui.plus_dt(0)

halt = True