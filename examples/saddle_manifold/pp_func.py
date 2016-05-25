from __future__ import division, absolute_import, print_function
# itertools, operator used for _filter_consecutive function
import itertools, operator
import os

from PyDSTool import *
from PyDSTool.errors import PyDSTool_ValueError
from PyDSTool.ModelContext import *
from PyDSTool.utils import findClosestPointIndex
from PyDSTool.common import args, metric, metric_L2, metric_weighted_L2, \
     metric_float, remain, fit_quadratic, fit_exponential, fit_diff_of_exp, \
     smooth_pts, nearest_2n_indices, make_poly_interpolated_curve, simple_bisection
from PyDSTool.common import _seq_types, _num_types
from PyDSTool.core.context_managers import RedirectStdout

import numpy as np
try:
    from numpy import unique
except ImportError:
    # older version of numpy
    from numpy import unique1d as unique
try:
    import matplotlib.pyplot as pp
except ImportError:
    pp = None

from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
from scipy.optimize import fsolve, minpack
from scipy.optimize import root, zeros
try:
    newton_meth = minpack.newton
except AttributeError:
    # newer version of scipy
    newton_meth = zeros.newton
from scipy import linspace, isfinite, sign, alltrue, sometrue, arctan, arctan2

from random import uniform
import copy
import sys
import six

norm = np.linalg.norm


from PyDSTool.Toolbox.phaseplane import *
from PyDSTool.Toolbox.phaseplane import bisection # not in __all__

# -----------------------
from fovea.diagnostics import diagnostic_manager
global dm
dm = diagnostic_manager('saddle_dm')
from fovea.graphics import gui
plotter = gui.plotter

def find_saddle_manifolds(fp, xname, ds=None, ds_gamma=None, ds_perp=None, tmax=None,
                          max_arclen=None, ic=None, eps=None, ev_dirn=1,
                          ic_ds=None, max_pts=1000, directions=(1,-1),
                          which=('s', 'u'), other_pts=None, rel_scale=None,
                          ds_perp_fac=0.75, verboselevel=0, fignum=None):
    """Compute any branch of the stable or unstable sub-manifolds of a saddle.
    Accepts fixed point instances of class fixedpoint_2D.

    Required inputs:
      fp:       fixed point object
      xname:    coordinate name of the x-axis variabe (e.g., for correct
                 orientation of verbose plotting)
      ds:       arc-length step size (**fixed**)
      ds_gamma: determines the positions of the Gamma_plus and Gamma_minus
                event surfaces (can be a real scalar or a pair if not symmetric)
      ds_perp:  initial perturbation from the local linear sub-manifolds to
                find starting points, computed in the direction of the eigenvalue.
      tmax:     maximum time to compute a test trajectory before 'failing' to find
                the Gamma event surface.
      max_arclen:  maximum arc length to compute
      max_pts:  maximum number of points to compute on each sub-manifold branch
      ic / ic_ds:
        Specify either ic or ic_ds for initial point (e.g. to restart the calc
        after a previous failure) or a certain distance from the saddle point.
      ev_dirn:  +1/-1. Event detection direction for Gamma_plus event.
          This may need to be flipped after trying one direction and getting errors
          that the event was not detected. An automated way to set this is not yet
          available, so you have to use trial and error or some forethought!
          Rule is: event direction code = 1 if, along the computed trajectory traj:
             gamm_ev(traj(ev_t - delta)) < 0 and gamma_ev(traj(ev_t + delta)) > 0
          for event detection time ev_t and small delta>0. The eigenvector along the
          flow towards the event surfaces determines which is "before" and which is
          "after" the event surface. (Also note that time will be reversed when
          computing the unstalbe manifold, which you have to take into account.)

    Optional inputs:
      eps: epsilon tolerance for manifold calculations (defaults to 1/100 times
          that of the FP objects passed to the function)

      which:  which sub-manifold to compute 's', 'u' or ('s', 'u').
        Default is both.

      directions:  which directions along chosen sub-manifolds? (1,), (-1,)
        or (1,-1). Default is both.

      rel_scale:  a pair giving relative scalings of x and y coordinates in
        the plane, to improve stepping in the different directions.
        e.g. (1,10) would make ds steps in the y-direction 10 times larger than
        in the x-direction. Default is (1,1)

      other_pts can be a list of points whose proximity will be checked,
    and the computation halted if they get within ds of the manifold.

      ds_perp_fac:  For advanced use only. If you get failures saying ds_perp
         too small and that initial displacement did not straddle manifold, try
         increasing this factor towards 1 (default 0.75). Especially for
         unstable manifolds, initial values for ds_perp may diverge, but if
         ds_perp is shrunk too quickly with this factor the sweet spot may be
         missed.

      verboselevel: 0 means silent, 1 means basic text info, 2 means extended info
         and also diagnostic plots. In diagnostic plots, the two solid black lines
         show the gamma plus/minus exit event lines; the blue trajectories are the
         test trajectories starting at the green cross test points;
         flow direction from IC shown as solid red line, its normal as dotted red.
      fignum:   Select figure number (defaults to 1)


     Returns:
       Dictionary keyed by 's' and 'u', containing dictionaries keyed by directions
         1 or -1 (ints) to a Pointset (if computed) or None (if not computed);
         parameterized by arc length, depending on user's selections for 'directions'
         argument and 'which' argument.
         E.g., if directions=(1,) and which = ('s',), returned structure looks like
            {'s': {1: <<pointset>>, -1: None}, 'u': {1: None, -1: None}}

    """
    if verboselevel > 1:
        figure_name, layer_name = plotter.active_layer
        _, layer_struct = plotter.active_layer_structs
        assert layer_struct is not None
    assert fp.classification == 'saddle' and not fp.degenerate
    if fp.evals[0] < 0:
        eval_s = fp.evals[0]
        eval_u = fp.evals[1]
        evec_s = fp.evecs[0]
        evec_u = fp.evecs[1]
    else:
        eval_s = fp.evals[1]
        eval_u = fp.evals[0]
        evec_s = fp.evecs[1]
        evec_u = fp.evecs[0]
    gen = fp.gen
    assert 'Gamma_out_plus' in gen.eventstruct, "Detection event surface(s) not present"
    assert 'Gamma_out_minus' in gen.eventstruct, "Detection event surface(s) not present"
    if eps is None:
        # Dividing fixed point's inherited epsilon tolerance by 100
        eps = fp.eps / 100
    ds_perp_eps = 1e-12
    if ds_perp_fac >= 1 or ds_perp_fac <= 0:
        raise ValueError("ds_perp_fac must be between 0 and 1")
    normord = fp.normord
    if rel_scale is None:
        rel_scale = (1,1)
    dsscaled = dx_scaled_2D(ds, rel_scale)
    if isinstance(ds_gamma, dict):
        assert len(ds_gamma) == 2, "Invalid value for ds_gamma"
        assert remain(list(ds_gamma.keys()), [1,-1]) == [], \
               "Invalid value for ds_gamma"
    else:
        try:
            ds_gamma = {1: ds_gamma, -1: ds_gamma}
        except:
            raise TypeError("Invalid type for ds_gamma")

    try:
        xcoord_ix = fp.point.coordnames.index(xname)
    except ValueError:
        raise ValueError("Invalid x coordinate name '%s'"%xname)
    else:
        # x coordinate index is either 0 or 1 for this 2D system
        # y coordinate index is therefore 1-xcoord_ix
        ycoord_ix = 1-xcoord_ix
    yname = fp.point.coordnames[ycoord_ix]

    if verboselevel>1:
        # validate coord names
        xn, yn = layer_struct.axes_vars
        if xname != xn and yname != yn:
            raise ValueError("x and y name mismatch with Plotter")

    def test_fn(x, dircode):
        if verboselevel>1:
            dm.log.msg("Integrate from test point", x=x[xname], y=x[yname], direction=dircode)
        gen.set(ics=x)
        try:
            test = gen.compute('test', dirn=dircode)
        except KeyboardInterrupt:
            raise
        except:
            raise RuntimeError("Integration failed")
        events = gen.getEvents()

        if verboselevel>1:
            pts=test.sample(coords=x.coordnames)
            # only show first 25 points unless Gamma bd not met
            plotter.add_data((pts[xname][:25],pts[yname][:25]), style='b-',
                            layer=layer_name,
                            name=dm.get_unique_name('test_traj_first25_'))

        if events['Gamma_out_plus'] is None:
            if events['Gamma_out_minus'] is None:
                if verboselevel>1:
                    pts = test.sample(coords=x.coordnames)
                    dm.log.msg("Error", err_msg="Did not reach Gamma surfaces",
                               status="fail", last_computed_point=pts[-1],
                               last_computed_time=pts['t'][-1])
                    plotter.add_data((pts[xname],pts[yname]), style='b-',
                                    layer=layer_name,
                                    name=dm.get_unique_name('test_traj_full'),
                                    log=dm.log)
                raise RuntimeError("Did not reach Gamma surfaces")
            else:
                # hit Gamma_out_minus
                if verboselevel>1:
                    dm.log.msg("Reached Gamma minus", t=events['Gamma_out_minus']['t'][0],
                               last_computed_point=pts[-1],
                               last_computed_time=pts['t'][-1])
                sgn = -1
        else:
            if events['Gamma_out_minus'] is None:
                # hit Gamma_out_plus
                if verboselevel>1:
                    dm.log.msg("Reached Gamma plus", t=events['Gamma_out_plus']['t'][0],
                               last_computed_point=pts[-1],
                               last_computed_time=pts['t'][-1])
                sgn = 1
            else:
                # both were non-None, i.e. both events happened: impossibru!
                if verboselevel>1:
                    pts = test.sample(coords=x.coordnames)
                    dm.log.msg("Error", err_msg="Both Gamma surfaces reached",
                               status="fail", last_computed_point=pts[-1],
                               last_computed_time=pts['t'][-1])
                    plotter.add_data((pts[xname],pts[yname]), style='b-',
                                    layer=layer_name,
                                    name=dm.get_unique_name('universe_fail'),
                                    log=dm.log)
                raise RuntimeError("Both Gamma surfaces reached, impossibly")
        return sgn

    def onto_manifold(x_ic, dn, normal_dir, dircode='f'):
        try:
            return bisection(test_fn, x_ic+dn*normal_dir, x_ic-dn*normal_dir,
                             args=(dircode,), xtol=eps, maxiter=100,
                             normord=normord)
        except AssertionError:
            if verboselevel>1:
                xp = x_ic+dn*normal_dir
                xm = x_ic-dn*normal_dir
                dm.log.msg("Error", err_msg="onto_manifold bisection fail",
                           status="fail", point_p=xp, point_m=xm)
                plotter.add_data([xp[xname],xp[yname]], style='gx',
                                 layer=layer_name,
                                 name=dm.get_unique_name('xp'), log=dm.log)
                plotter.add_data([xm[xname],xm[yname]], style='gx',
                                 layer=layer_name,
                                 name=dm.get_unique_name('xm'), log=dm.log)
                plotter.show()
            raise RuntimeError("ds_perp too small? +/- initial displacement did not straddle manifold")
        except RuntimeError:
            if verboselevel>1:
                xp = x_ic+dn*normal_dir
                xm = x_ic-dn*normal_dir
                dm.log.msg("Error", err_msg="onto_manifold bisection fail",
                           status="fail", point_p=xp, point_m=xm)
                plotter.add_data([xp[xname],xp[yname]], style='gx',
                                 layer=layer_struct.name,
                                 name=dm.get_unique_name('xp'), log=dm.log)
                plotter.add_data([xm[xname],xm[yname]], style='gx',
                                 layer=layer_struct.name,
                                 name=dm.get_unique_name('xm'), log=dm.log)
                plotter.show()
            raise

    gen.eventstruct['Gamma_out_plus'].activeFlag=True  # terminal
    gen.eventstruct['Gamma_out_minus'].activeFlag=True  # terminal
    assert tmax > 0
    manifold = {'s': {1: None, -1: None}, 'u': {1: None, -1: None}}
    man_names = {'s': 'stable', 'u': 'unstable'}

    for w in which:
        # w = 's' => stable branch
        # w = 'u' => unstable branch
        if verboselevel>0:
            print("Starting %s branch" % man_names[w])
        if w == 's':
            col = 'g'
            w_sgn = -1
            integ_dircode = 'f'
            evec = evec_u
            evec_other = evec_s
        elif w == 'u':
            col = 'r'
            w_sgn = 1
            integ_dircode = 'b'
            evec = evec_s
            evec_other = evec_u
        # set Gamma_out surfaces on "outgoing" branch
        # (polarity is arbitrary)
        p0_plus = fp.point + ds_gamma[1]*evec
        p0_minus = fp.point - ds_gamma[-1]*evec
        evec_perp = get_perp(evec)
        gen.eventstruct.setEventDir('Gamma_out_plus', ev_dirn)
        gen.eventstruct.setEventDir('Gamma_out_minus', -ev_dirn)
        gen.set(pars={'Gamma_out_plus_p_'+xname: p0_plus[xname],
                      'Gamma_out_plus_p_'+yname: p0_plus[yname],
                      'Gamma_out_plus_dp_'+xname: evec_perp[xname],
                      'Gamma_out_plus_dp_'+yname: evec_perp[yname],
                      'Gamma_out_minus_p_'+xname: p0_minus[xname],
                      'Gamma_out_minus_p_'+yname: p0_minus[yname],
                      'Gamma_out_minus_dp_'+xname: evec_perp[xname],
                      'Gamma_out_minus_dp_'+yname: evec_perp[yname],
                      ##                  'fp_'+xname: fp.point[xname], 'fp_'+yname: fp.point[yname]
                      },
                tdata = [0,tmax])
        if verboselevel>1:
            if fignum is None:
                fignum=figure()
            else:
                figure(fignum)
            # plot event surfaces for gamma plus and minus exit events
            # ISSUE: Convert to plotter.add_data
            plot([p0_plus[xname]-dsscaled*evec_perp[xname],p0_plus[xname]+dsscaled*evec_perp[xname]],
                 [p0_plus[yname]-dsscaled*evec_perp[yname],p0_plus[yname]+dsscaled*evec_perp[yname]], 'k-', linewidth=2)
            plot([p0_minus[xname]-dsscaled*evec_perp[xname],p0_minus[xname]+dsscaled*evec_perp[xname]],
                 [p0_minus[yname]-dsscaled*evec_perp[yname],p0_minus[yname]+dsscaled*evec_perp[yname]], 'k-', linewidth=2)
            draw()
        check_other_pts = other_pts is not None
        if ic_ds is None:
            ic_ds = dsscaled
        else:
            ic_ds = dx_scaled_2D(ic_ds, rel_scale)
        if ic is None:
            ic = fp.point
            f_ic = -w_sgn * evec_other
            dirn_fix = 1 # not used for this case
            if verboselevel>0:
                # ISSUE: Convert to log entry
                print("f_ic from evec_other")
                print("evec_other " + str(evec_other))
                print("f_ic = " + str(f_ic))
            curve_len = 0
            # initial estimate x0 = a point close to f.p. along manifold with
            # opposite stability
        else:
            # initial curve length from previous independent variable, if present
            # otherwise, assume zero
            if isinstance(ic, Pointset):
                assert len(ic) == 1, "Only pass a length-1 pointset"
                # (guarantee curve_len > 0)
                # BUG: for direction=-1 case, arc_len will be negative
                # and index 0 will have the smallest arc_len, not the
                # largest. Better not to use ic as Pointset option and
                # fix arc_len outside of call
                curve_len = abs(ic['arc_len'][0])
                ic = ic[0]
            else:
                curve_len = 0
            # ensure correct sign relative to starting point (if ic is None)
            sgns_orig = sign(-w_sgn * evec_other)
            f_ic_alpha = gen.Rhs(0, ic, gen.pars)  # array in alpha order
            # f_ic here isn't normalized to length 1 like the case above that uses
            # evec_other (which is already normalized)
            f_ic = Point({xname: f_ic_alpha[xcoord_ix], yname: f_ic_alpha[ycoord_ix]})
            sgns_f_ic = sign(f_ic)
            if any(sgns_orig != sgns_f_ic):
                dirn_fix = -1
                f_ic = -f_ic
            else:
                dirn_fix = 1
            if verboselevel>0:
                # ISSUE: Convert to log entry
                print("f_ic = " + str(f_ic))
        for sgn in directions:
            piece = {}
            if verboselevel>0:
                # ISSUE: Convert to log entry
                print("Starting direction", sgn)
            # PREDICTION
            x0_ic = ic+w_sgn*sgn*ic_ds*f_ic/norm(f_ic, normord)
            if verboselevel>1:
                figure(fignum)
                # show starting point (initial estimate) as green circle
                # ISSUE: Convert to plotter.add_data
                plot(x0_ic[xname], x0_ic[yname], 'go', linewidth=1)
            # put x0 initial estimate onto stable manifold
            f_alpha = dirn_fix * gen.Rhs(0, x0_ic, gen.pars)  # array in alpha order
            f = Point({xname: f_alpha[xcoord_ix], yname: f_alpha[ycoord_ix]})
            normf = norm(f, normord)
            norm_to_flow = get_perp(f/normf)
            if verboselevel>1:
                # show flow direction from IC as solid red line
                plotter.add_data(([x0_ic[xname], x0_ic[xname]+dsscaled*f[xname]/normf],
                     [x0_ic[yname], x0_ic[yname]+dsscaled*f[yname]/normf]),
                     style='r-', name=dm.get_unique_name('flow_fwd'), log=dm.log)
                # show normal to flow direction from IC as dotted red line
                plotter.add_data(([x0_ic[xname], x0_ic[xname]+dsscaled*norm_to_flow[xname]],
                     [x0_ic[yname], x0_ic[yname]+dsscaled*norm_to_flow[yname]]),
                     style='r:', name=dm.get_unique_name('flow_perp'), log=dm.log)
            ds_perp_default = ds_perp
            # CORRECTION
            while ds_perp > ds_perp_eps:
                try:
                    x = onto_manifold(x0_ic, ds_perp, norm_to_flow,
                                      dircode=integ_dircode)
                except RuntimeError as e:
                    ds_perp *= ds_perp_fac
                else:
                    break
            if ds_perp <= ds_perp_eps:
                # RuntimeError was raised and could not continue reducing ds_perp
                print("ds_perp reached lower tolerance =", ds_perp_eps)
                print(e)
                raise RuntimeError("Initial point did not converge")
            else:
                curve_len += norm(x-ic, normord)
                piece[sgn*curve_len] = x
                num_pts = 1
                last_x = x
                if verboselevel>0:
                    print("Initial point converged to (%.6f, %.6f)\n" % \
                          (x[xname], x[yname]))
            ds_perp = ds_perp_default
            last_f = f_ic
            # step backwards along local linear flow to predict next starting
            # position on manifold
            while curve_len < max_arclen and num_pts < max_pts:
                if verboselevel>0:
                    # ISSUE: Convert to plotter.add_data
                    figure(fignum)
                    plot(last_x[xname], last_x[yname], col+'.', linewidth=1)
                if check_other_pts and sometrue([norm(last_x - pt, normord) < ds \
                                                 for pt in other_pts]):
                    # we've hit a different fixed point (or other feature), so stop
                    break
                f_alpha = dirn_fix * gen.Rhs(0, last_x, gen.pars) # array
                f = Point({xname: f_alpha[xcoord_ix], yname: f_alpha[ycoord_ix]})
                if all(sign(f) != sign(last_f)):
                    f = -f
                    # on other side of manifold so must keep stepping in the
                    # same direction, therefore switch signs!
                # PREDICTION
                x_ic = last_x + w_sgn*sgn*dsscaled*f/norm(f,normord)
                last_f = f
                if verboselevel>1:
                    print("\nStarting from point ", last_x)
                    delta = w_sgn*sgn*dsscaled*f/norm(f,normord)
                    print("Trying point ", x_ic, "in direction (%.6f, %.6f)\n" % (delta[xname], delta[yname]))
                ds_perp = ds_perp_default
                # CORRECTION
                while ds_perp > ds_perp_eps:
                    try:
                        x = onto_manifold(x_ic, ds_perp, get_perp(f/norm(f,normord)),
                                          dircode=integ_dircode)
                    except RuntimeError as e:
                        ds_perp *= 0.75
                    else:
                        break
                if ds_perp <= ds_perp_eps:
                    # RuntimeError was raised and could not continue reducing ds_perp
                    print("ds_perp reached lower tolerance =", ds_perp_eps)
                    print(e)
                    break  # end while search
                else:
                    curve_len += norm(x-last_x, normord)
                    piece[sgn*curve_len] = x
                    last_x = x
                    num_pts += 1
                    if verboselevel>1:
                        print("\nManifold has %i points" % num_pts)
                    elif verboselevel>0:
                        print(".", end=' ')
                        sys.stdout.flush()
            indepvar, piece_sorted = sortedDictLists(piece, byvalue=False)
            manifold[w][sgn] = pointsToPointset(piece_sorted, indepvarname='arc_len',
                                                indepvararray=indepvar, norm=normord)
        if verboselevel>0:
            # finish the line on stdout
            print(" ")
    gen.eventstruct['Gamma_out_plus'].activeFlag=False
    gen.eventstruct['Gamma_out_minus'].activeFlag=False
##    gen.eventstruct['fp_closest'].activeFlag=False

    return manifold


