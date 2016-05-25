"""
"""
from __future__ import division
import PyDSTool as dst
from PyDSTool.Toolbox import phaseplane as pp
import numpy as np
from matplotlib import pyplot as plt

import fovea
import fovea.graphics as gx
from fovea.graphics import tracker
from fovea.graphics import gui

import pp_func
global dm
# reuse existing objects
dm = pp_func.dm
plotter = gui.plotter
dm.use_dir('saddle_logs')
dm.make_log('log.json')

# connect the diagnostic_manager to the plotter
plotter.dm = dm

import tinydb as tdb

# user selection
print("Trying to use Radau integrator, but use Vode")
print("  if you don't have the C integrators working")
gentype = dst.Generator.Radau_ODEsystem

all_plots = False

# ------------------------------

def build_sys():
    # we must give a name
    DSargs = dst.args(name='M345_A3_Bead_on_a_rotating_hoop')

    # parameters
    DSargs.pars = {'g': 0,
                   'd': 0.3}

    # rhs of the differential equation
    DSargs.varspecs = {'phi': 'nu',
                       'nu': '-d*nu + g*sin(phi)*cos(phi) - sin(phi)'}

    # initial conditions
    DSargs.ics = {'phi': 0, 'nu': 0}

    # set the domain of integration.
    # (increased domain size to explore around phi=-pi saddle)
    DSargs.xdomain = {'phi': [-2*np.pi, 2*np.pi], 'nu': [-4, 4]}

    # allow tdomain to be infinite, set default tdata here
    DSargs.tdata = [0, 50]

    # to avoid typos / bugs, use built-in Symbolic differentation!
    f = [DSargs.varspecs['phi'], DSargs.varspecs['nu']]
    Df=dst.Diff(f, ['phi', 'nu'])
    DSargs.fnspecs = {'Jacobian': (['t','phi','nu'],
                                   str(Df.renderForCode()))}
    # yields """[[0, 1], [g*cos(phi)*cos(phi) - g*sin(phi)*sin(phi) - cos(phi), -d]]""")}
    print("Jacobian computed as:\n" + str(Df.renderForCode()))

    # Make auxiliary functions to define event lines near saddle
    res = pp.make_distance_to_line_auxfn('Gamma_out_plus',
                                      'Gamma_out_plus_fn',
                                      ('phi','nu'), True)
    man_pars = res['pars']
    man_auxfns = res['auxfn']
    res = pp.make_distance_to_line_auxfn('Gamma_out_minus',
                                      'Gamma_out_minus_fn',
                                      ('phi','nu'), True)
    man_pars.extend(res['pars'])
    man_auxfns.update(res['auxfn'])

    # update param values with defaults (0)
    for p in man_pars:
        DSargs.pars[p] = 0

    if gentype in [dst.Generator.Vode_ODEsystem, dst.Generator.Euler_ODEsystem]:
        targetlang = 'python'
    else:
        targetlang = 'c'

    DSargs.fnspecs.update(man_auxfns)
    ev_plus = dst.Events.makeZeroCrossEvent(expr='Gamma_out_plus_fn(%s,%s)'%('phi','nu'),
                                         dircode=0,
                                         argDict={'name': 'Gamma_out_plus',
                                                  'eventtol': 1e-5,
                                                  'eventdelay': 1e-3,
                                                  'starttime': 0,
                                                  'precise': False,
                                                  'active': False,
                                                  'term': True},
                                         targetlang=targetlang,
                                         varnames=['phi','nu'],
                                         fnspecs=man_auxfns,
                                         parnames=man_pars
                                        )
    ev_minus = dst.Events.makeZeroCrossEvent(expr='Gamma_out_minus_fn(%s,%s)'%('phi','nu'),
                                         dircode=0,
                                         argDict={'name': 'Gamma_out_minus',
                                                  'eventtol': 1e-5,
                                                  'eventdelay': 1e-3,
                                                  'starttime': 0,
                                                  'precise': False,
                                                  'active': False,
                                                  'term': True},
                                         targetlang=targetlang,
                                         varnames=['phi','nu'],
                                         fnspecs=man_auxfns,
                                         parnames=man_pars
                                         )

    DSargs.events = [ev_plus, ev_minus]

    # an instance of the 'Generator' class.
    print("Initializing generator...")
    return gentype(DSargs)


def plot_PP_fp(fp, plotter_layer, coords=None, do_evecs=True, markersize=10):
    """Draw fixed point from phaseplane.fixedpoint_2D objects.

    Optional do_evecs (default True) draws eigenvectors
    """
    figure, layer = plotter.active_layer
    fig_struct, layer_struct = plotter.active_layer_structs

    if coords is None:
        x, y = layer_struct.axes_vars #fp.fp_coords # not stored in used order, just alphabetic order
    else:
        x, y = coords
    # When matplotlib implements half-full circle markers
    #if fp.classification == 'saddle':
        # half-open circle
    if fp.stability == 'u':
        style = 'wo'
    elif fp.stability == 'c':
        style = 'co'
    else: # 's'
        style = 'ko'

    # Add point and eigenvectors
    name = 'fp'
    dm.log.msg('Plot FP', name=name, layer=plotter_layer)
    plotter.add_point(pp.Point2D(saddle.point[x], saddle.point[y]),
                     style=style, layer=plotter_layer, name=name)

    if do_evecs:
        for i, evec in enumerate(fp.evecs):
            fp_pt = fp.point
            evalue = fp.evals[i]
            if fp.classification == 'saddle':
                # this labeling will only work for a saddle
                if np.real(evalue) < 0:
                    name = 'evec_s'
                    style = 'g-.'
                else:
                    name = 'evec_u'
                    style = 'r-.'
            else:
                if np.real(evalue) < 0:
                    name = 'evec_s'
                    style = 'g-.'
                else:
                    name = 'evec_u'
                    style = 'r-.'
            p1, p2 = fovea.graphics.force_line_to_extent(fp_pt, fp_pt+pp.Point2D(evec, xname=x, yname=y),
                                                         dict(zip(layer_struct.axes_vars,layer_struct.scale)),
                                                         layer_struct.axes_vars)
            line = np.array((p1, p2))
            name = dm.get_unique_name(name)
            dm.log.msg('Plot evec', name=name, layer=plotter_layer)
            plotter.add_data(line.T, style=style, layer=plotter_layer, name=name)


def plot_PP_vf_custom(gen, xname, yname, N=20, subdomain=None, scale_exp=0):
    """Draw 2D vector field in (xname, yname) coordinates of given Generator,
    sampling on a uniform grid of n by n points.

    Optional subdomain dictionary specifies axes limits in each variable,
    otherwise Generator's xdomain attribute will be used.

    For systems of dimension > 2, the non-phase plane variables will be held
      constant at their initial condition values set in the Generator.

    Optional scale_exp is an exponent (domain is all reals) which rescales
      size of arrows in case of disparate scales in the vector field. Larger
      values of scale magnify the arrow sizes. For stiff vector fields, values
      from -3 to 3 may be necessary to resolve arrows in certain regions.

    Requires matplotlib 0.99 or later
    """
    assert N > 1
    xdom = gen.xdomain[xname]
    ydom = gen.xdomain[yname]
    if subdomain is not None:
        try:
            xdom = subdomain[xname]
        except KeyError:
            pass
        try:
            ydom = subdomain[yname]
        except KeyError:
            pass
    assert all(dst.isfinite(xdom)), "Must specify a finite domain for x direction"
    assert all(dst.isfinite(ydom)), "Must specify a finite domain for y direction"
    w = xdom[1]-xdom[0]
    h = ydom[1]-ydom[0]

    xdict = gen.initialconditions.copy()

    xix = gen.funcspec.vars.index(xname)
    yix = gen.funcspec.vars.index(yname)

    xs = np.linspace(xdom[0], xdom[1], N)
    ys = np.linspace(ydom[0], ydom[1], N)

    X, Y = np.meshgrid(xs, ys)
    dxs, dys = np.meshgrid(xs, ys)

    dz_big = 0
    vec_dict = {}

    for xi, x in enumerate(xs):
        for yi, y in enumerate(ys):
            xdict.update({xname: x, yname: y})
            dx, dy = gen.Rhs(0, xdict)[[xix, yix]]
            # note order of indices
            dxs[yi,xi] = dx
            dys[yi,xi] = dy
            dz = np.linalg.norm((dx,dy))
            if dz > dz_big:
                dz_big = dz

    plt.quiver(X, Y, dxs, dys, angles='xy', pivot='middle', units='inches',
               scale=dz_big*max(h,w)/(10*np.exp(2*scale_exp)), lw=0.01/np.exp(scale_exp-1),
               headwidth=max(2,1.5/(np.exp(scale_exp-1))),
               #headlength=2*max(2,1.5/(exp(scale_exp-1))),
               width=0.001*max(h,w), minshaft=2, minlength=0.001)

    ax = plt.gca()

    print("xdom: ", xdom)
    print("ydom: ", ydom)
    ax.set_xlim(xdom)
    ax.set_ylim(ydom)
    plt.draw()


def test_traj(ic, tend=20):
    """Convenience function for exploring trajectories"""
    ode_sys.set(ics=ic, tdata=[0,tend])
    traj = ode_sys.compute('test')
    pts = traj.sample()
    plt.plot(pts['phi'], pts['nu'], 'k:', lw=1)
    plt.plot(pts['phi'][0], pts['nu'][0], 'ko')




# -----------------------------------

ode_sys = build_sys()

# domain of interest
DOI = ([-3.25,-3],[-0.075,0.075])
plotter.clean()
plotter.add_fig('Master',
               title='Manifold Computation Diagnostic',
               xlabel='phi', ylabel='nu', # redundant vs. arrange_fig options below?
               domain=DOI)

common_layer_setup = {'axes_vars': ['phi', 'nu'],
                      'scale': DOI}

plotter.add_layer('fp_data', **common_layer_setup)
plotter.add_layer('vf_data', **common_layer_setup)
plotter.add_layer('nullcline_data', **common_layer_setup)
plotter.add_layer('manifold_data', **common_layer_setup)
plotter.add_layer('manifold_metadata', **common_layer_setup)

# phase plane tools are in the Toolbox module aliased as 'pp'

# plot vector field, using a scale exponent to ensure arrows are well spaced
# and sized
if all_plots:
    plot_PP_vf_custom(ode_sys, 'vf_data', scale_exp=-0.25)

plotter.arrange_fig([1,1], {'11': {'name': 'PP',
                                  'scale': DOI,
                                  'layers': ['fp_data', 'manifold_data', 'manifold_metadata'],
                                  'axes_vars': ['phi', 'nu']}
                           })

# find fixed points
fp_coords = pp.find_fixedpoints(ode_sys, eps=1e-6)

# n=3 uses three starting points in the domain to find nullcline parts, to an
# accuracy of eps=1e-8, and a maximum step for the solver of 0.1 units.
# The fixed points found is also provided to help locate the nullclines.
if all_plots:
    nulls_x, nulls_y = pp.find_nullclines(ode_sys, 'phi', 'nu', n=3,
                                      eps=1e-6, max_step=0.1, fps=fp_coords)

# plot the nullclines
if all_plots:
    plt.plot(nulls_x[:,0], nulls_x[:,1], 'b')
    plt.plot(nulls_y[:,0], nulls_y[:,1], 'g')

# plot the fixed points
fps = []
for fp_coord in fp_coords:
    fps.append( pp.fixedpoint_2D(ode_sys, dst.Point(fp_coord)) )

saddle = fps[1]

plotter.set_active_layer('fp_data')
plot_PP_fp(saddle, 'fp_data', do_evecs=True, markersize=7)

gui.build_plotter((8,8), with_times=False)

# magBound change ensures quicker determination of divergence during
# manifold computations. max_pts must be larger when we are further
# away from the fixed point.
ode_sys.set(algparams={'magBound': 10000})

def plot_manifold(man, which, style='k.-'):
    for sgn in (-1, 1):
        if man[which][sgn] is not None:
            print("There were %i points in sub-manifold %s, direction %i" % (len(man[which][sgn]), which, sgn))
            plt.plot(man[which][sgn]['phi'], man[which][sgn]['nu'], style)

# KEY to manifold parts dictionary key codes
# u means unstable branch
# s means stable branch
# +1 means 'forwards' direction
# -1 means 'backwards' direction

manifold_parts = {
            'u': {1: dst.Pointset(indepvarname='arc_len',
                              indepvararray=[0],
                              coorddict=saddle.point),
                  -1: dst.Pointset(indepvarname='arc_len',
                              indepvararray=[0],
                              coorddict=saddle.point)},
            's': {1: dst.Pointset(indepvarname='arc_len',
                                          indepvararray=[0],
                                          coorddict=saddle.point),
                 -1: dst.Pointset(indepvarname='arc_len',
                                          indepvararray=[0],
                                          coorddict=saddle.point)}
            }

ode_sys.set(ics=saddle.point, tdata=[0,60])

which_man = 's'
dirn = -1
verbose = 2
man_part = manifold_parts[which_man][dirn]

ds_perp = 0.02 # initial value for stage 2
ds_gamma = 0.3 # initial value for stage 2

# for speed, and because of symmetry, only compute long arcs on one side
if dirn == 1:
    max_arclen = 4
else:
    max_arclen = 9
max_pts = 1


plotter.set_active_layer('manifold_metadata')


# first stage (only called once)
print("  First stage...")

print("Note: f_ic does not seem to correspond to eigenvector direction")
ode_sys.set(algparams={'max_pts': 20000})
def test(ic_ds, ds_perp, ds_gamma, eps, ev_dirn, tmax):
    # for first, single step, ds is irrelevant: only ic_ds is important
    man_new = pp_func.find_saddle_manifolds(saddle, 'phi', ds=0.1, ds_gamma=ds_gamma,
                                   ds_perp=ds_perp, tmax=tmax, max_arclen=max_arclen, eps=eps,
                                   ic_ds=ic_ds, max_pts=1, directions=(dirn,), ev_dirn=ev_dirn,
                                   which=(which_man,), other_pts=[fps[0].point, fps[2].point],
                                   rel_scale=(1,1), verboselevel=verbose, fignum=1)
    return man_new

man_new = test(0.01, 0.001, 0.01, 1e-3, 1, 10)

print(dm.db.search(tdb.where('event')=="Added plot data"))

plotter.show()
print(plotter.shown)
plt.show()

##man_new = pp.find_saddle_manifolds(saddle, 'phi', ds=0.004, ds_gamma=0.02,
##                ds_perp=0.005, tmax=60, max_arclen=max_arclen, eps=2e-5,
##                ic_ds=0.0002, max_pts=max_pts, directions=(dirn,), ev_dirn=1,
##                which=(which_man,), other_pts=[fps[0].point, fps[2].point],
##                rel_scale=(1,1), verboselevel=verbose, fignum=1)
part = man_new[which_man][dirn]
if dirn == 1:
    select = -1
else:
    # backwards means that points will be added behind the initial point
    # in terms of Pointset order, so the "last" point computed will be
    # at index 0
    select = 0
part.indepvararray += man_part.indepvararray[select]
man_part.insert(part)

# stop here
1/0

# Stage two (called repeatedly until length while loop satisfied)
print("  Continuing in stage 2...")
ode_sys.set(algparams={'max_pts': 100000})
if dirn == 1:
    select = -1
else:
    # backwards means that points will be added behind the initial point
    # in terms of Pointset order, so the "last" point computed will be
    # at index 0
    select = 0
attempt_num += 1
try:
    # groups of max_pts/4 at a time
    man_new = pp.find_saddle_manifolds(saddle, 'phi', ds=0.06, ds_gamma=ds_gamma,
                ds_perp=ds_perp, tmax=40, max_arclen=max_arclen, eps=2e-5,
                ic=man_part[select], max_pts=int(max_pts/4.0),
                directions=(dirn,), ev_dirn=1,
                which=(which_man,), other_pts=[fps[0].point, fps[2].point],
                rel_scale=(1,1), verboselevel=verbose, fignum=1)
except RuntimeError:
    # proceed with what we've got
    print("Initial convergence error: Proceeding with what we've got!")
    ds_perp *= 2
    ds_gamma *= 4
else:
    part = man_new[which_man][dirn]
    part.indepvararray += man_part.indepvararray[select]
    man_part.insert(part)

plot_manifold(manifold_parts, 's', 'k-')
plot_manifold(manifold_parts, 'u', 'r-')

plt.show()
