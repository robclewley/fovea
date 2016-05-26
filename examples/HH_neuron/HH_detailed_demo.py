"""
This is the main run script for the detailed demo involving
Hodgkin-Huxley analysis.
"""

from PyDSTool.Toolbox.dssrt import *
from PyDSTool.Toolbox.phaseplane import *
import PyDSTool as dst
import numpy as np
import scipy as sp

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

##def test_fp(V):
##    pt = {'V': V,
##          'K.n': gen.auxfns.K_dssrt_fn_ninf(V)+0.001,
##          'Na.m': gen.auxfns.Na_dssrt_fn_minf(V)}
##    omega_n = da_rate_signed.calc_psi('K.n', pt)
##    omega_m = da_rate_signed.calc_psi('Na.m', pt)
##    # Vdot will be zero
##    return omega_n + omega_m

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
        # gui, gen and da_reg_signed are globals
        self.loc = local_linear_2D(gen, da_reg_signed, 'V', xvar, gui.points[gui.ix],
                      with_exit_auxfn=False, tol=1e-5)


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
        loc = self.loc
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
        plotter.add_data([[pt['Na.m'], pt['Na.m']+dm_dt*self.vel_arrow_scale],
                         [pt['V'], pt['V']+dV_dt*self.vel_arrow_scale]],
                layer='state_vel_mV', name='state', style=vel_vec_style, force=True)

        dvinf_dt_n = ptsDSSRT[gui.ix]['Omega_n']
        # self.dQ_dt('vinf', gui.ix, gui.points)
        plotter.add_data([[pt['Na.m'], pt['Na.m']],
                         [pt['vinf'], pt['vinf']+dvinf_dt_n*self.vel_arrow_scale]],
                layer='state_vel_mV', name='vinf', style=vinf_vec_style, force=True)

        if self.first_call:
            plotter.add_data([gui.points['Na.m'], gui.points['V']],
                            layer='vfp_mV', name='traj', style='y')

        # Virtual fixed point and linearized nullclines
        if 'fast_m' in model.name:
            vfp = None
            with_jac = False
            do_fps = False
            fast_vars = ['Na.m']
        else:
            loc.analyze(pt) # this is a slow step!
            vfp = loc.fp
            with_jac = False #True
            do_fps = False #True
            fast_vars = None
            lin_ms = np.linspace(sc[0][0], sc[0][1], 3)
            lin_vinfs = [loc.lin.auxfns.vinf(m) for m in lin_ms]
            lin_minfs = [loc.lin.auxfns.fnim(m) for m in lin_ms]
            plotter.add_data([lin_ms, lin_vinfs], layer='vfp_mV', name='lin_nullV', style='b:', force=True)
            plotter.add_data([lin_ms, lin_minfs], layer='vfp_mV', name='lin_nullm', style='g:', force=True)

        # update (or create) points
        try:
            plotter.set_point('state_pt', Point2D(pt['Na.m'], pt['V']), 'points_mV')
            plotter.set_point('vinf_pt', Point2D(pt['Na.m'], pt['vinf']), 'points_mV')
            if vfp:
                plotter.set_point('vfp_pt', Point2D(vfp['m'], vfp['v']), 'points_mV')
        except KeyError:
            plotter.add_point(Point2D(pt['Na.m'], pt['V']),
                         layer='points_mV', style='ko', name='state_pt')
            plotter.add_point(Point2D(pt['Na.m'], pt['vinf']),
                         layer='points_mV', style='bx', name='vinf_pt')
            if vfp:
                plotter.add_point(Point2D(vfp['m'], vfp['v']), layer='points_mV',
                             name='vfp_pt', style={'color': 'y', 'marker': 'o',
                                                   'markersize': 5})

        d = fig_struct.layers['nullclines_mV'].data

        if not preComputed and gui._mouseUp:
            print("\nComputing phase plane...")
            print("  Current time = %.4f" % (time))

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
            plotter.add_data(self.nully, layer='nullclines_mV', style=self.nullcY_style,
                            name='yNull_'+str(time), force=force)

            dssrt_data = ptsDSSRT[gui.ix]
            # ensure don't plot "time left to horizon collision" line if distance is negative
            if dssrt_data['d'] > 0:
                state = (pt['Na.m'], pt['V'])
                d_state = (dssrt_data['d']*sin(dssrt_data['c']), dssrt_data['d']*cos(dssrt_data['c']))
                plotter.add_data([[state[0]+d_state[0], state[0]],
                              [state[1]+d_state[1], state[1]]], layer='horiz_PP',
                            style=horiz_style, name='horiz_'+str(time), force=force)

            # delete update 'wait' notice
            ax.texts = []
            #ax.clear()
            gui.clear_axes(ax)

            if only_var is None:
                # nullx is added second so will be the second line
                self.nullx = castNullArray(nulls['nullcX'])
                plotter.add_data(self.nullx, layer='nullclines_mV',
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
            print("  Phase plane rebuild completed.\n")
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
        loc = self.loc
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
        plotter.add_data([[pt['K.n'], pt['K.n']+dn_dt*self.vel_arrow_scale],
                         [pt['V'], pt['V']+dV_dt*self.vel_arrow_scale]],
                layer='state_vel_nV', name='state', style=vel_vec_style, force=True)

        dvinf_dt_m = ptsDSSRT[gui.ix]['Omega_m']
        # self.dQ_dt('vinf', gui.ix, gui.points)
        plotter.add_data([[pt['K.n'], pt['K.n']],
                         [pt['vinf'], pt['vinf']+dvinf_dt_m*self.vel_arrow_scale]],
                layer='state_vel_nV', name='vinf', style=vinf_vec_style, force=True)

        if self.first_call:
            plotter.add_data([gui.points['K.n'], gui.points['V']],
                            layer='vfp_nV', name='traj', style='y')
            plotter.add_data([gui.points['K.n'], gui.points['vinf']],
                            layer='vfp_nV', name='quasiVnull', style='m--')

        vs = np.linspace(sc[1][0], sc[1][1], 100)
        x = dict(pt).copy()

        def vinf(n, v):
            x['K.n'] = n
            x['V'] = v
            x['Na.m'] = gen.auxfns.Na_dssrt_fn_minf(v)
            # assume autonomous system
            return model.Rhs(0, x, asarray=False)['V']

        vinfs_inv_n = np.array([fsolve(vinf, gen.auxfns.K_dssrt_fn_ninf(v), args=(v,)) for v in vs]).T[0]
        if self.first_call:
            plotter.add_data([vinfs_inv_n, vs], layer='vfp_nV', name='vinf_fastm', style='b--')
        else:
            plotter.set_data('vfp_nV', data={'vinf_fastm': {'data': [vinfs_inv_n, vs], 'style':'b--', 'display': True}}, display=True)

        # Virtual fixed point and linearized nullclines
        if 'fast_m' in model.name:
            vfp = None
            with_jac = False
            do_fps = False
            fast_vars = ['Na.m']
        else:
            loc.analyze(pt)
            vfp = loc.fp
            with_jac = False #True
            do_fps = False #True
            fast_vars = None
            lin_ns = np.linspace(sc[0][0], sc[0][1], 3)
            lin_vinfs = [loc.lin.auxfns.vinf(n) for n in lin_ns]
            lin_ninfs = [loc.lin.auxfns.fnin(n) for n in lin_ns]
            plotter.add_data([lin_ns, lin_vinfs], layer='vfp_nV',
                            name='lin_nullV', style='b:', force=True)
            plotter.add_data([lin_ns, lin_ninfs], layer='vfp_nV',
                            name='lin_nulln', style='r:', force=True)

        # update (or create) points
        try:
            plotter.set_point('state_pt', Point2D(pt['K.n'], pt['V']), 'points_nV')
            plotter.set_point('vinf_pt', Point2D(pt['K.n'], pt['vinf']), 'points_nV')
            if vfp:
                plotter.set_point('vfp_pt', Point2D(vfp['n'], vfp['v']), 'points_nV')
        except KeyError:
            plotter.add_point(Point2D(pt['K.n'], pt['V']),
                         layer='points_nV', style='ko', name='state_pt')
            plotter.add_point(Point2D(pt['K.n'], pt['vinf']),
                         layer='points_nV', style='bx', name='vinf_pt')
            if vfp:
                plotter.add_point(Point2D(vfp['n'], vfp['v']), layer='points_nV',
                             name='vfp_pt', style={'color': 'y', 'marker': 'o',
                                                   'markersize': 5})

        d = fig_struct.layers['nullclines_nV'].data

        if not preComputed and gui._mouseUp:
            print("\nComputing phase plane...")
            print("  Current time = %.4f" % (time))

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
            plotter.add_data(self.nully, layer='nullclines_nV', style=self.nullcY_style,
                            name='yNull_'+str(time), force=force)

            # delete update 'wait' notice
            ax.texts = []
            #ax.clear()
            gui.clear_axes(ax)

            if only_var is None:
                # nullx is added second so will be the second line
                self.nullx = castNullArray(nulls['nullcX'])
                plotter.add_data(self.nullx, layer='nullclines_nV',
                                style=self.nullcX_style,
                                name='xNull', force=force)

            #if force:
            #    rescale = sc
            #else:
            #    rescale = None
            gui.build_layers(['nullclines_nV', 'points_nV', 'state_vel_nV', 'vfp_nV'],
                            ax, rescale=sc, figure='Master')

            self.last_scale = sc
            print("  Phase plane rebuild completed.\n")
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


# rough prototype functionality to add "feature"-like data to existing plots
def add_epochs_to_layer(kind, to_layer, epochs, style, figure=None):
    """
    Add DSSRT epoch times to target layer as cross points in specified color.
    Target layer may not be a 'dynamic' layer.

    kind          string e.g. 'Psi' or 'Omega'
    to_layer      e.g., 'A'
    """
    fig_struct, figure = plotter._resolve_fig(figure)
    ep_times = [ep.t0 for ep in epochs]

    if to_layer not in fig_struct.layers:
        raise ValueError("Target layer %s not found" % to_layer)
    layer_name = kind + ' epochs @ ' + to_layer
    plotter.add_layer(layer_name, kind='epochs_'+kind)
    layer = fig_struct.layers[to_layer]
    for traj_name, traj in layer.trajs.items():
        if layer.kind == 'data':
            vals = traj(ep_times)['y']
            plotter.add_data([ep_times, vals], layer=layer_name,
                            style=style,
                            name=traj_name+'.epochs')
    return layer_name



def animate_frames(t0, t1, n_frames=100,
                   fname_base=None,
                   format='png'):
    if fname_base is None:
        'screenshots_gK%i/fovea_gK%i_' % (Kgmax, Kgmax)

    assert t0 < t1
    for i, t in enumerate(np.linspace(t0, t1, n_frames)):
        print(i)
        sys.stdout.flush()
        gui.set_time(t)
        plt.savefig(fname_base+str(i+1)+'.'+format, format=format)

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

pert_time = 4.5  # choose > t_end for no pert
dV = 0.2

##test_ic = {'K.n': 0.37220277852490802,
##           'Na.m': 0.080387043479386036,
##           'V': -59.5}
##model.set(ics=test_ic)
##force_DSSRT = True
force_DSSRT = False #True

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
    if pert_time < t_end:
        model.set(tdata=[0, pert_time])
        model.compute('ref1')
        ref_traj1 = model['ref1']
        pts = ref_traj1.sample()
        new_ic = ref_traj1(pert_time)
        new_ic['V'] += dV
        model.set(ics=new_ic,
                  tdata=[0, t_end-pert_time])
        model.compute('ref2')
        ref_traj2 = model['ref2']
        pts2 = ref_traj2.sample()
        pts2.indepvararray += pert_time
        pts.extend(pts2, skipMatchingIndepvar=True)
        ref_traj = dst.pointset_to_traj(pts)
        model.set(ics=orig_ics)
    else:
        model.set(tdata=[0, t_end])
        model.compute('ref')
        ref_traj = model['ref']
else:
    # get periodic orbit
    t_end = 20
    ref_traj, ref_pts, ref_tmin, ref_tmax = do_traj(model, t_end,
                                                    do_plot=False)

# re-sample traj at constant dt and declare to GUI
trajPts = ref_traj.sample(dt=0.01)[:-40] #  cheap way to avoid overlap from pts not being periodic
#gui.add_data_points(trajPts)


## ----- ----- ----- ----- ----- ----- ##
## CREATE DIAGNOSTIC OBJECT            ##
## ----- ----- ----- ----- ----- ----- ##

plotter.clean()
plotter.add_fig('Master', title='Geometric Dynamic Analysis: '+dssrt_name,
               tdom=[0, t_end])

# Add layers and their data

plotter.add_layer('V')
plotter.add_data([trajPts['t'], trajPts['V']], layer='V', style='k-',
                name='V')
plotter.add_data([trajPts['t'], trajPts['vinf']], layer='V', style='k:',
                name='Vinf')

plotter.add_layer('activs')
plotter.add_data([trajPts['t'], trajPts['Na.m']], layer='activs', style='g-',
                name='m')
plotter.add_data([trajPts['t'], trajPts['Na.minf']], layer='activs', style='g--',
                name='minf')

plotter.add_data([trajPts['t'], trajPts['K.n']], layer='activs', style='r-',
                name='n')
plotter.add_data([trajPts['t'], trajPts['K.ninf']], layer='activs', style='r--',
                name='ninf')

plotter.add_data([trajPts['t'], trajPts['tauv']], layer='activs', style='b:',
                name='tauv')
plotter.add_data([trajPts['t'], trajPts['Na.taum']], layer='activs', style='g:',
                name='taum')
plotter.add_data([trajPts['t'], trajPts['K.taun']], layer='activs', style='r:',
                name='taun')

## ----- ----- ----- ----- ----- ----- ##
## CALCULATE DSSRT (OMEGAs, PSIs)      ##
## ----- ----- ----- ----- ----- ----- ##

(ptsDSSRT, strDSSRT, epochs_reg, epochs_rate), (da_reg, da_rate, \
                                    da_reg_signed, da_rate_signed) = \
    getHH_DSSRT(model, trajPts, dssrt_name, header=header, verbose=3,
                force=force_DSSRT)
# strDSSRT is a table for display purposes and can be saved to .csv if needed

# print out epochs to screen
print("\nEpochs regular (psi):")
show_epochs(epochs_reg)

print("\nEpochs rate (omega):")
show_epochs(epochs_rate)


##### Build layers, sub-plots
ptsDSSRT['Psi_n'][0] = 1 # prevents NaN messing up plot
ptsDSSRT['Omega_n'][0] = 1 # prevents NaN messing up plot

plotter.add_layer('A')
plotter.add_data([ptsDSSRT['t'], ptsDSSRT['A']], layer='A', style='k-',
                name='A')
plotter.add_data([ptsDSSRT['t'], ptsDSSRT['accn']], layer='A', style='k:',
                name='accn')
plotter.add_layer('Omegas')
plotter.add_data([ptsDSSRT['t'], ptsDSSRT['Omega_n']], layer='Omegas', style='r-',
                name='Omega_n')
plotter.add_data([ptsDSSRT['t'], ptsDSSRT['Omega_m']], layer='Omegas', style='g-',
                name='Omega_m')
plotter.add_layer('Vdot')
plotter.add_data([ptsDSSRT['t'], ptsDSSRT['Vdot']], layer='Vdot', style='b-',
                name='Vdot')

plotter.add_layer('A_hline')
plotter.add_hline(0, layer='A_hline', style='k:', name='zero')

plotter.add_layer('rel times')
plotter.add_data([ptsDSSRT['t'], ptsDSSRT['tleft']], layer='rel times',
                style='g-', name='tleft')
plotter.add_data([ptsDSSRT['t'], ptsDSSRT['horiz_t']], layer='rel times',
                style='b-', name='horiz_t')
plotter.add_layer('t_hline')
plotter.add_hline(0, layer='t_hline', style='k:', name='zero')


##plotter.add_layer('dom_m_to_n')
##plotter.add_data([ptsDSSRT['t'], ptsDSSRT['Psi_m']/ptsDSSRT['Psi_n']],
##                layer='dom_m_to_n', style='k-', name='psi_rat')
##plotter.add_data([ptsDSSRT['t'], ptsDSSRT['Omega_m']/ptsDSSRT['Omega_n']],
##                layer='dom_m_to_n', style='b-', name='omega_rat')
##plotter.add_layer('rat_hlines')
##plotter.add_hline(1, layer='rat_hlines', style='k:', name='plus_one')
##plotter.add_hline(-1, layer='rat_hlines', style='k:', name='minus_one')


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
    plotter.add_layer(PP_layer_name, dynamic=True)
    if xvar == 'Na.m':
        plotter.add_layer('horiz_PP')
        nullcX_style = 'g-'
        PPclass = PPcallback_m
    else:
        # no horizon layer for K.n
        nullcX_style = 'r-'
        PPclass = PPcallback_n
    PPplot = PPclass(xvar, nullcY_style = {'color': 'b', 'linestyle': '-', 'linewidth': 1},
                                                      nullcX_style=nullcX_style)
    gui.dynamicPlotFns[PP_layer_name] = PPplot
    plotter.add_layer('points_'+suffix)
    plotter.add_layer('state_vel_'+suffix)
    plotter.add_layer('vfp_'+suffix)

make_layer('Na.m')
make_layer('K.n')

# sub-plot specs: (row, col) integer coords start at top left
dPlot11 = {'name': 'Trajectory',
           'layers': ['V'],
           'scale': [None, [-72, 40]],
           'axes_vars': ['t', 'V'] }

dPlot12 = {'name': 'Activations, Time scales',
           'layers': ['activs'],
           'scale': [None, [0,1]],
           'axes_vars': ['t', 'no units, ms'] }

dPlot23 = {'name': 'Magic A, Omegas',
           'layers': ['A', 'A_hline', 'Omegas', 'Vdot'],
           'scale': [None, [-2,2]],
           'axes_vars': ['t', 'mV / ms'] }

dPlot13 = dict(name='Linearized event times',
           layers=['rel times', 't_hline'],
           axes_vars=['t', 'time'],
           scale=[None, [-2,20]])

pp1_name = 'Na-V Phaseplane'
pp1_dom = [[0,0.25], [-75,40]] #-35]]
pp1_vars = ['m', 'V']
pp2_name = 'K-V Phaseplane'
# focus on first knee
##if 'typeI' in model.name and 'typeII' not in model.name:
##    pp2_dom = [[0.,0.3], [-75,-35]]
##else:
##    pp2_dom = [[0.32,0.5], [-75,-35]]
pp2_dom = [[0.2,1], [-75,60]]
pp2_vars = ['n', 'V']

dPlot21 = {'name': pp1_name,
           'scale': pp1_dom,
           'layers': ['nullclines_mV', 'horiz_PP', 'points_mV', 'state_vel_mV', 'vfp_mV'],
           'axes_vars': pp1_vars}

dPlot22 = {'name': pp2_name,
           'scale': pp2_dom,
           'layers': ['nullclines_nV', 'points_nV', 'state_vel_nV', 'vfp_nV'],
           'axes_vars': pp2_vars}

##dPlot22 = {'name': 'Dominance m:n',
##           'layers': ['dom_m_to_n', 'rat_hlines'],
##           'scale': [None, [-5, 5]],
##           'axes_vars': ['t', 'Psi ratio']}

dPlot_dict = {'11': dPlot11, '12': dPlot12, '21': dPlot21, '22': dPlot22,
                           '13': dPlot13, '23': dPlot23}

# add epoch info to first declared layer in each time sub-plot
for ixstr, dP in dPlot_dict.items():
    if dP['axes_vars'][0] != 't':
        continue
    # could improve by checking which layers are lines and have kind='data'
    data_layer = dP['layers'][0]
    epoch_layer_name_psi = add_epochs_to_layer('Psi', data_layer, epochs_reg, 'kx')
    epoch_layer_name_omega = add_epochs_to_layer('Omega', data_layer, epochs_rate, 'b+')
    dP['layers'].extend([epoch_layer_name_psi, epoch_layer_name_omega])

plotter.arrange_fig([2,3], dPlot_dict)

gui.build_plotter((14,7))


##import os
##if 'WINGDB_ACTIVE' in os.environ:
##    plt.show()

plotter.show_legends(subplot='Times')

##plt.figure(2)
##plt.plot(ptsDSSRT['t'], ptsDSSRT['theta'])
##plt.plot(ptsDSSRT['t'], ptsDSSRT['a'])

plt.show()
