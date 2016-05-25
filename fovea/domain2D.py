"""
Domain2D tools

Depends on the Shapely and descartes packages

"""
from __future__ import division, absolute_import

import PyDSTool.Toolbox.phaseplane as pp
import PyDSTool as dst

import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import shapely.geometry as geom
import warnings
from shapely.geometry import polygon as P
from descartes.patch import PolygonPatch
import numpy as np

class polygon_domain(object):
    """
    """
    def __init__(self, c, p1, condition_func, rtol=0.1,
                 nsides=10, edge_len_rtol=2, max_radius_fac=500,
                 reseed_rate_fac=10):
        """
        c is the center seed point
        p1 is any point that defines the initial radius
           relative to c.
           (p1 will become one vertex of the polygon.)
        condition_func is a function of (x,y) tuple that
           returns a scalar, or a sequence of such functions.
           The domain is defined by the set of connected points
           with the same sign as condition_func(c).
        rtol (default 0.1) is the relative tolerance of
           the growth stopping condition, as a fraction of the
           initial "radius" of the regular polygon defined
           by |p1-c|.
        nsides (default 10) is the initial number of sides
           for the polygons.
        edge_len_rtol (default 2) is the factor of |p1-c| that
           determines the maximum side length of a polygon edge
           before new vertex is added adaptively during domain
           growth.
        max_radius_fac (default 1000) is the factor of |p1-c|
           that determines the largest distance from c that
           a grown polygon can be before the iterations stop.
        reseed_rate_fac (default 10) is the factor of |p1-c|
           that determines how far from the original center
           seed point the polygon can grow in any direction
           before a new seed polygon is added ???????????????
        """
        self.initial_center = np.asarray(c)
        self.initial_p = np.asarray(p1)
        self.initial_r = self.initial_p - self.initial_center
        if callable(condition_func):
            # single function
            target_fsign = np.sign(condition_func(c))
            self.boolfunc = lambda p: np.sign(condition_func(p, fsign = target_fsign)) == \
                                         target_fsign
        else:
            # list of functions
            target_fsigns = [np.sign(f(c)) for f in condition_func]
            self.boolfunc = lambda p: [np.sign(f(p)) for f in condition_func] == \
                                         target_fsigns
        # check that p1 satifies same sign
        if not self.boolfunc(p1):
            #raise dst.PyDSTool_ValueError("f(p1) has different sign to f(c)")\
            #print("Error: f(p1) has different sign to f(c)")
            warnings.warn("Warning: f(p1) has different sign to f(c)")
        else:
            # distance scale and initial trust radius
            self.r0 = np.linalg.norm(p1-c)
            self.max_radius = max_radius_fac*self.r0
            self.edge_resolution_atol = rtol*self.r0
            self.edge_length_atol = edge_len_rtol*self.r0
            # angle of trust radius line to x-axis
            theta1 = math.acos(self.initial_r[0]/self.r0)
            dtheta = 2*math.pi/nsides
            # distribute angles uniformly over the circle
            # and compute nsides-1 new polygon vertices from p1
            pt_list = [self.initial_center + \
                               self.r0*np.array([math.cos(theta1+dtheta*i),
                                                 math.sin(theta1+dtheta*i)])
                               for i in range(nsides)]
            self.polygon = P.Polygon(pt_list)
            self.nsides = nsides

            self.stop_growing = [False]*self.nsides

    def grow(self, max_iters=300):
        for i in range(max_iters):
            self.grow_step()

    def grow_step(self, verbose=False):
        # first and last pt are identical
        pts = np.array(self.polygon.exterior.coords)[:-1]
        c = self.initial_center
        dr = 0.1*self.r0
        f = self.boolfunc
        new_pts = pts.copy()  # default stay the same
        for i, p in enumerate(pts):
            if self.stop_growing[i]:
                continue
            r_unit = (p-c)/np.linalg.norm(p-c)
            # test new point
            pnew = p + dr*r_unit
            if np.linalg.norm(pnew-c) > self.max_radius:
                self.stop_growing = [True]*self.nsides
                if verbose:
                    print("All points stopped growing because too big")
                break
            elif f(pnew):
                # accept point
                new_pts[i] = pnew
            else:
                # bisect old p with pnew
                # for now, define success here
                self.stop_growing[i] = True
                if verbose:
                    print("Point %i stopped growing because success"%i)

        self.polygon = P.Polygon(new_pts)

        # check that new side lengths are not too long
        edge_len_diffs = np.asarray((edge_lengths(self.polygon) - self.edge_length_atol) > 0, dtype='int')
        bad_edge_ixs = np.argwhere(edge_len_diffs).T.flatten()
        if len(bad_edge_ixs) > 0:
            final_pts = new_pts.copy()
            for ix in bad_edge_ixs:
                # Vertices for edge in polygon object are ix and ix+1 which is guaranteed
                # to be in range of indices because last point is repeated.

                # IDEAL CASE:
                # Scan along radius 1/2 way between the two vertices to find a safe point
                # to start a new polygon with same trust radius as this polygon started
                # with

                # PRACTICAL CASE:
                # Revert these indices to their previous values and stop growing them --
                # relying on user to start a different polygon and merge
                final_pts[ix] = pts[ix]
                ix2 = np.mod(ix+1,self.nsides)
                final_pts[ix2] = pts[ix2]
                self.stop_growing[ix] = self.stop_growing[ix2] = True
                if verbose:
                    print("Points %i and %i stopped growing because edge too long" % (ix, ix2))
                    print("   ", pts[ix], pts[ix2])
            self.polygon = P.Polygon(final_pts)


def edge_lengths(polyg):
    pts = np.array(polyg.exterior.coords)
    # first and last points are identical so ignore last point
    return np.array([np.linalg.norm(dp) for dp in (pts[1:-1] - pts[:-2])])

def merge_to_polygon(p1, p2):
    """
    Merges two polygon objects, two polygon_domain objects,
    or combinations thereof, and returns a polygon object
    consisting of the merger of the polygons of the two
    argument objects
    """
    if hasattr(p1, 'polygon'):
        # assume polygon_domain
        p1 = p1.polygon
    if hasattr(p2, 'polygon'):
        # assume polygon_domain
        p2 = p2.polygon
    return p1.union(p2)

class GUI_domain_handler(object):
    def __init__(self, gui, verbose=False):
        self.func = None
        self.domain = None  # polygon object
        self.gui = gui  # hook back (was for Bombardier GUIrocket object)
        self.polygon_domain_obj = None
        # grow states per polygon:
        #  0 = not begun, no domain object created
        #  1 = center placed, no initial radius
        #  2 = initial polygon defined, domain object created
        #  3 = final polygon defined
        #  4 = 3+1 (new center)
        #  5 = 3+2
        #  (6 = new polygon defined) --> 3 after merge (not yet implemented)
        self.gui_grow_state = 0
        self.center_pt = None
        self.p1_pt = None
        self.verbose = verbose

    def assign_criterion_func(self, func):
        """
        func( (x,y) ) --> real scalar

        func is a zero-crossing function whose nullcline
        defines the boundary of a domain, or a sequence of
        such functions
        """
        self.func = func

    def event(self, evcode):
        """
        Simple finite state machine
        legal event codes:
          'key' = domain key pressed
          'mouse' = mouse primary button pressed
          'reset' = clear data and revert to state 0
                    (keep any existing domain and function)
        returns success/validity of state update from event
        """
        if evcode == 'reset':
            self.gui_grow_state = 0
            if self.verbose:
                print("State 0 reset")
            self.polygon_domain_obj = None
            self.unshow_domain()
            return True
        elif evcode == 'key':
            if self.gui_grow_state > 0:
                # key doesn't make sense here, but GUIrocket.key_on
                # logic should prevent arriving here
                if self.verbose:
                    print("Domain grow state already active: ignoring keypress")
                return False
            else:
                # state 0
                if self.verbose:
                    print("State 0: setting up 1")
                self.gui.mouse_cid = self.gui.fig.canvas.mpl_connect('button_release_event', self.mouse_event_make_dom_c)
                return True
        elif evcode == 'mouse':
            if self.gui_grow_state == 0:
                if self.verbose:
                    print("State 0 -> 1")
                self.gui_grow_state = 1
                return True
            elif self.gui_grow_state == 1:
                if self.verbose:
                    print("State 1 -> 2")
                self.gui_grow_state = 2
                return True
            else:
                if self.verbose:
                    print("Not yet supported")
                return False

    def mouse_event_make_dom_c(self, ev):
        if ev.inaxes not in self.gui.cb_axes:
            print('Must select axes for which callbacks have been defined.')
            return

        self.curr_axes = ev.inaxes
        if self.verbose:
            print("In make_dom_c")
        # release mouse event control
        self.gui.fig.canvas.mpl_disconnect(self.gui.mouse_cid)
        # update state
        self.event('mouse')
        if self.gui_grow_state != 1:
            if self.verbose:
                print("make_dom_c failed")
            self.gui.mouse_wait_state_owner = None
            return
        # assign to c
        self.center_pt = pp.Point2D(ev.xdata, ev.ydata)
        # display c
        #self.gui.selected_object_temphandle = self.gui.ax.plot(ev.xdata, ev.ydata, 'go')[0]
        self.gui.selected_object_temphandle = ev.inaxes.plot(ev.xdata, ev.ydata, 'go')[0]
        self.gui.fig.canvas.draw()
        # switch control to make_dom_p1
        self.gui.mouse_cid = self.gui.fig.canvas.mpl_connect('button_release_event', self.mouse_event_make_dom_p1)

    def mouse_event_make_dom_p1(self, ev):
        if self.verbose:
            print("In make_dom_p1")
        # release mouse event control
        self.gui.fig.canvas.mpl_disconnect(self.gui.mouse_cid)
        self.gui.mouse_wait_state_owner = None
        # update state
        self.event('mouse')
        if self.gui_grow_state != 2:
            if self.verbose:
                print("make_dom_p1 failed")
            return
        # assign to p1
        self.p1_pt = pp.Point2D(ev.xdata, ev.ydata)
        # delete display of c
        self.gui.selected_object_temphandle.remove()
        # create then iterate polygon_domain object
        self.polygon_domain_obj = polygon_domain(self.center_pt,
                                                     self.p1_pt,
                                                     self.func,
                                                     nsides=40,
                                                     edge_len_rtol=3)
        if self.verbose:
            print("Growing domain")
        try:
            self.polygon_domain_obj.grow()
            self.domain = self.polygon_domain_obj.polygon
            self.gui_grow_state = 3
            self.show_domain(ev)
            if self.verbose:
                print("Domain complete")
        except AttributeError:
            return

    def show_domain(self, ev):
        fig_struct, figure = self.gui.plotter._resolve_fig(None)

        ##ISSUE: This try/except should probably be moved to domain_GUI.
        try:
            self.gui.plotter._resolve_layer(figure, 'gx_objects')
        except KeyError:
            self.gui.plotter.add_layer('gx_objects', subplot=ev.inaxes, kind = 'obj')
            print("Created layer gx_objects to support domain polygon")

        xs, ys = self.polygon_domain_obj.polygon.exterior.xy
        self.gui.set_selected_object(fovea.graphics.domain_GUI(self.gui, [xs, ys], name= 'domain1', subplot= ev.inaxes), figure= self.plotter.currFig)
        self.gui.plotter.show(ignore_wait = True)

    def unshow_domain(self):
        if self.gui.selected_object_temphandle is not None:

            try:
                self.gui.selected_object_temphandle.remove()
            except ValueError:
                # sequence
                try:
                    for th in self.gui.selected_object_temphandle:
                        th.remove()
                except TypeError:
                    pass

        self.gui.fig.canvas.draw()


#Local import
#Must import at bottom, as graphics already imports domain2D.
import fovea.graphics