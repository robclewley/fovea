"""
Rocket multi-body simulation inspired by
Bombardiers' Guild mobile app game
"""

from __future__ import division
import os
from PyDSTool import *
import PyDSTool.Toolbox.phaseplane as pp
import PyDSTool as dst   # for potentially generalizable functions and classes to use
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RectangleSelector

import fovea
import fovea.domain2D as dom
from fovea import common, prep
import fovea.graphics as gx


gentype = 'vode'

# let rocket mass be negligible
# and G chosen to absorb m_rocket

# grav constant
G = 35

# Generic scale separation threshold
scale_thresh = 3.0

# graphics / control params
xdomain_halfwidth = .7  # should derive from YAML setup in calling script
maxangle = 80 # degrees
da_dict = dict(zip( ('h','H','j','J','n','m'),
                    (-1, -10, 1, 10, -0.1, 0.1)))
dv_dict = dict(zip( ('d','D','f','F','c','v'),
                    (-1, -10, 1, 10, -0.1, 0.1)))

print("Change angle keys:")
print(da_dict)

print("Change velocity keys:")
print(dv_dict)

# other keys used in GUIrocket:
# l = make a line of interest (click-drag-release)
# SPACE = measure forces at clicked mouse point
# s = snap clicked mouse point to closest point on trajectory
# . (o) = grow a 2D domain
# g = GO! (run simulation)

dom_key = '.'
change_mouse_state_keys = ['l', 's', ' '] + [dom_key]


# Disable mpl interactive key defaults
plt.rcParams['keymap.save'] = 'ctrl+s' # s
plt.rcParams['keymap.xscale'] = '' # l
plt.rcParams['keymap.yscale'] = '' # L
plt.rcParams['keymap.back'] = ['left', 'backspace'] # c
plt.rcParams['keymap.forward'] = 'right' # v
plt.rcParams['keymap.zoom'] = '' # o
plt.rcParams['keymap.grid'] = ''  # g

# initial value
next_fighandle = 1

_non_pickle_attr = ['ax', 'fig',
                    'trajline', 'startpt', 'endpt', 'quartiles',
                    'widgets',
                    'RS_line']
# attributes that themselves contain non-picklable objects
_non_pickle_subattr = ['context_objects',  # e.g. line_GUI
                       'tracked_objects']  # e.g. measure with dynamic fns


class GUIrocket(gx.diagnosticGUI):
    def __init__(self, bodies, title, axisbgcol='black'):
        """
        bodies is a dict-like mapping of 1 or more:
           <ID int>: {'density': <float>,
                      'radius': <float>,
                      'position': [<float>, <float>]}
        """
        global next_fighandle

        plotter = gx.Plotter()
        gx.diagnosticGUI.__init__(self, plotter)

        self.current_domain_handler = dom.GUI_domain_handler(self)

        # --- SPECIFIC TO BOMBARDIER
        # Setup shoot params
        self.vel = 0.8
        self.ang = 0
        self.da = 0.005
        self.dv = 0.0005
        # used for color-coding trajectories by speed
        self.maxspeed = 2.2

        # one time graphics setup
        # for plot handles
        self.trajline = None
        self.startpt = None
        self.endpt = None
        self.quartiles = None

        # Axes background colour
        self.axisbgcol = axisbgcol

        #Setup code
        DOI = [(-xdomain_halfwidth,xdomain_halfwidth),(0,1)]
        self.clean() # in case rerun in same session
        self.add_fig('master',
                        title='Bombardier',
                        xlabel='x', ylabel='y',
                        domain=DOI)

        #Setup all layers

        self.add_layer('trajs')
        self.add_layer('bodies', kind='patch')
        self.add_layer('text', kind='text')

        self.name = 'gamespace'

        self.setup({'11':
                    {'name': self.name,
                     'scale': DOI,
                     'layers':['trajs', 'bodies', 'text'],
                     'callbacks':'*',
                     'axes_vars': ['x', 'y']
                     }
                    },
                  size=(9, 7), with_times=False, basic_widgets=False)

        self.fignum = 1

        fig_struct, fig = self.plotter._resolve_fig('master')
        self.ax = fig_struct.arrange['11']['axes_obj']

        self.add_widget(Slider, callback=self.updateAng, axlims = (0.1, 0.055, 0.65, 0.03),
                      label='Shoot Angle', valmin= -maxangle, valmax= maxangle,
                      valinit= self.ang, color='b', dragging=False, valfmt='%2.3f')

        self.add_widget(Slider, callback=self.updateVel, axlims=(0.1, 0.02, 0.65, 0.03),
                      label='Shoot Speed', valmin=0.01, valmax=2,
                      valinit=self.vel, color='b',
                      dragging=False, valfmt='%1.4f')


        # assume max of N-2 planetoid bodies + target + source
        self.N = len(bodies)
        self.gen_versioner = common.gen_versioner(os.path.abspath('.'),
                                                         self.name,
                                                         'simgen_N%i'%self.N,
                                                         gentype, 1)

        # Make this more generic for ABC
        self.setup_pars(bodies)

        # --- END OF BOMBARDIER SPECIFICS

        # Move these to a _recreate method than can be reused for un-pickling

        self.add_widget(Button, callback=self.go, axlims=(0.005, 0.1, 0.045, 0.03), label='Go!')

        # context_changed flag set when new objects created using declare_in_context(),
        # and unset when Generator is created with the new context code included
        self.context_changed = False
        self.setup_gen(self.model_namer)

        self.mouse_cid = None # event-connection ID
        self.go(run=False)
        # force call to graphics_refresh because run=False above
        self.graphics_refresh(cla=False)
        # TEMP
        #plt.show()

        # next_fighandle for whenever a new model is put in a new figure (new game instance)
        next_fighandle += 1


    def graphics_refresh(self, cla=True):
        if cla:
            self.ax.cla()

        #Make quartiles
        xquarts = Point({'x': 4})
        yquarts = Point({'y': 4})

        try:
            n = len(self.points)
            coorddict = {'xq':
                         {'x':'xq', 'y':'yq', 'layer':'trajs', 'name':'quarts1', 'style':'kd'}
                         }
            quarts = Pointset({'coordarray': np.array([[self.points['x'][int(0.25*n)], self.points['x'][int(0.5*n)], self.points['x'][int(0.75*n)]],
                                           [self.points['y'][int(0.25*n)], self.points['y'][int(0.5*n)], self.points['y'][int(0.75*n)]]]),
                      'coordnames': ['xq', 'yq']})
            self.add_data_points(quarts, coorddict=coorddict)

        except TypeError:
            pass

        #Traj Pointset
        coorddict = {'x':
                     {'x':'x', 'y':'y','layer':'trajs','name':'data1', 'object':'collection'},
                      #{'x':'x', 'y':'y','layer':'trajs', 'object':'collection'},
                     'speed':
                     {'map_color_to':'x'}
                     }
        self.add_data_points(self.points, coorddict=coorddict)

        #Bodies Pointset
        bodsPoints = Pointset({'coordarray': np.array([[self.pos[i][0] for i in range(len(self.pos))],
                                       [self.pos[i][1] for i in range(len(self.pos))],
                                       [self.radii[i] for i in range(len(self.radii))]]),
                  'coordnames': ['px', 'py', 'radii']})
        coorddict = {'px':
                     {'x':'px', 'y':'py','layer':'bodies','name':'bods1', 'style':'g', 'object':'circle'},
                     'radii':
                     {'map_radius_to':'px'}
                     }
        self.add_data_points(bodsPoints, coorddict=coorddict)

        pos = np.array(self.pos).transpose()
        for i in range(len(pos[0])):
            self.plotter.add_text(pos[0][i], pos[1][i], i, style='k', layer='text')

        self.plotter.show(rebuild=False)

    # Methods for pickling protocol
    def __getstate__(self):
        d = copy(self.__dict__)
        for fname, finfo in self._funcreg.items():
            try:
                del d[fname]
            except KeyError:
                pass
        # delete MPL objects
        for obj in some_list:
            try:
                del d['']
            except KeyError:
                pass
        return d

    def __setstate__(self, state):
        # INCOMPLETE!
        self.__dict__.update(state)
        self._stuff = None
        if something != {}:
            self._recreate() # or re-call __init__

    def _recreate(self):
        raise NotImplementedError

    #def declare_in_context(self, con_obj):
        ## context_changed flag set when new objects created and unset when Generator is
        ## created with the new context code included
        #self.context_changed = True
        #self.context_objects.append(con_obj)

    def __str__(self):
        return self.name

    def setup_pars(self, data):
        # Should generalize to non-bombardier application
        N = self.N
        radii = {}
        density = {}
        pos = {}
        for i, body in data.items():
            pos[i] = pp.Point2D(body['position'][0], body['position'][1],
                                labels={'body': i})
            radii[i] = body['radius']
            density[i] = body['density']
        ixs = range(N)
        self.radii = [radii[i] for i in ixs]
        self.density = [density[i] for i in ixs]
        self.pos = [pos[i] for i in ixs] # planet positions
        self.masses = [density[i]*np.pi*r*r for (i,r) in enumerate(self.radii)]
        rdict = dict([('r%i' %i, self.radii[i]) for i in ixs])
        mdict = dict([('m%i' %i, self.masses[i]) for i in ixs])
        posxdict = dict([('bx%i' %i, pos[i][0]) for i in ixs])
        posydict = dict([('by%i' %i, pos[i][1]) for i in ixs])
        pardict = {'G': G}  # global param for gravitational constant
        pardict.update(rdict)
        pardict.update(mdict)
        pardict.update(posxdict)
        pardict.update(posydict)
        self.body_pars = pardict
        self.icpos = np.array((0.0, 0.08))
        self.icvel = np.array((0.0, 0.0))

    def make_gen(self, pardict, name):
        # scrape GUI diagnostic object extras for generator
        extra_events = []
        extra_fnspecs = {}
        extra_pars = {}
        extra_auxvars = {}
        for con_obj in self.context_objects.values():
            extra_events.append(con_obj.extra_events)
            extra_fnspecs.update(con_obj.extra_fnspecs)
            extra_pars.update(con_obj.extra_pars)
            extra_auxvars.update(con_obj.extra_auxvars)

        Fx_str = ""
        Fy_str = ""
        for i in range(self.N):
            Fx_str += "-G*m%i*(x-bx%i)/pow(d(x,y,bx%i,by%i),3)" % (i,i,i,i)
            Fy_str += "-G*m%i*(y-by%i)/pow(d(x,y,bx%i,by%i),3)" % (i,i,i,i)

        DSargs = args()
        DSargs.varspecs = {'vx': Fx_str, 'x': 'vx',
                           'vy': Fy_str, 'y': 'vy',
                           'Fx_out': 'Fx(x,y)', 'Fy_out': 'Fy(x,y)',
                           'speed': 'sqrt(vx*vx+vy*vy)',
                           'bearing': '90-180*atan2(vy,vx)/pi'}
        DSargs.varspecs.update(extra_auxvars)
        auxfndict = {'Fx': (['x', 'y'], Fx_str),
                     'Fy': (['x', 'y'], Fy_str),
                     'd': (['xx', 'yy', 'x1', 'y1'], "sqrt((xx-x1)*(xx-x1)+(yy-y1)*(yy-y1))")
                    }
        DSargs.auxvars = ['Fx_out', 'Fy_out', 'speed', 'bearing'] + \
            list(extra_auxvars.keys())
        DSargs.pars = pardict
        DSargs.pars.update(extra_pars)
        DSargs.fnspecs = auxfndict
        DSargs.fnspecs.update(extra_fnspecs)
        DSargs.algparams = {'init_step':0.001,
                            'max_step': 0.01,
                            'max_pts': 20000,
                            'maxevtpts': 2,
                            'refine': 5}

        targetlang = \
            self.gen_versioner._targetlangs[self.gen_versioner.gen_type]

        # Events for external boundaries (left, right, top, bottom)
        Lev = Events.makeZeroCrossEvent('x+%f'%xdomain_halfwidth, -1,
                                        {'name': 'Lev',
                                         'eventtol': 1e-5,
                                         'precise': True,
                                         'term': True},
                                        varnames=['x'],
                                        targetlang=targetlang)
        Rev = Events.makeZeroCrossEvent('x-%f'%xdomain_halfwidth, 1,
                                        {'name': 'Rev',
                                         'eventtol': 1e-5,
                                         'precise': True,
                                         'term': True},
                                        varnames=['x'],
                                        targetlang=targetlang)
        Tev = Events.makeZeroCrossEvent('y-1', 1,
                                        {'name': 'Tev',
                                         'eventtol': 1e-5,
                                         'precise': True,
                                         'term': True},
                                        varnames=['y'],
                                        targetlang=targetlang)
        Bev = Events.makeZeroCrossEvent('y', -1,
                                        {'name': 'Bev',
                                         'eventtol': 1e-5,
                                         'precise': True,
                                         'term': True},
                                        varnames=['y'],
                                        targetlang=targetlang)

        # Events for planetoids
        bevs = []
        for i in range(self.N):
            bev = Events.makeZeroCrossEvent('d(x,y,bx%i,by%i)-r%i' % (i,i,i),
                                            -1,
                                        {'name': 'b%iev' %i,
                                         'eventtol': 1e-5,
                                         'precise': True,
                                         'term': True},
                                        varnames=['x','y'],
                                        parnames=list(pardict.keys()),
                                        fnspecs=auxfndict,
                                        targetlang=targetlang)
            bevs.append(bev)

        DSargs.events = [Lev, Rev, Tev, Bev] + bevs + extra_events
        DSargs.checklevel = 2
        DSargs.ics = {'x': self.icpos[0], 'y': self.icpos[1],
                      'vx': 0., 'vy': 1.5}
        DSargs.name = name
        DSargs.tdomain = [0, 10000]
        DSargs.tdata = [0, 50]

        # turns arguments into Generator then embed into Model object
        self.model = self.gen_versioner.make(DSargs)


    def go(self, run=True):
        """
        Note: This method can only start a trajectory from the
        launcher at the bottom of the screen!

        To shoot from a specific point that's been set by hand,
        call self.run() then self.graphics_refresh(cla=False)
        """
        a = self.ang
        v = self.vel
        # Angle a of shooting is relative to vertical, up to +/- maxangle degrees
        if a > maxangle:
            # assume is vestigial from a different initial condition
            a = maxangle
        elif a < -maxangle:
            a = -maxangle
        rad = pi*(a-90)/180.
        x = self.radii[0]*cos(rad)
        y = -self.radii[0]*sin(rad)
        vx = v*cos(rad)
        vy = -v*sin(rad)
        self.model.set(ics={'vx': vx, 'vy': vy,
                             'x': x, 'y': y})
        if run:
            self.run()
            self.graphics_refresh(cla=False)
        self.masterWin.canvas.draw()
        plt.draw()

    def set(self, pair, ic=None, by_vel=False):
        """Set solution pair (ang, speed) and optional (x,y)
        initial condition, where ang is in degrees.

        With option by_vel=True (default False),
         the pair will be treated as (vx, vy) instead
        """
        assert len(pair) == 2
        if ic is not None:
            assert 'x' in ic and 'y' in ic and len(ic) == 2
            self.model.set(ics=ic)
            self.icpos = ic
            if by_vel:
                vx, vy = pair
                # both conversions in this section are -90?
                self.ang = 180*atan2(vy,vx)/pi - 90
                self.vel = sqrt(vx*vx+vy*vy)
            else:
                # can't set ang and vel according to rules for regular
                # shooting because we are reconstructing a partial
                # trajectory out in space
                self.ang, self.vel = pair
                rad = pi*(self.ang-90)/180.
                vx = self.vel*cos(rad)
                vy = -self.vel*sin(rad)
            self.model.set(ics={'vx': vx, 'vy': vy})
        else:
            self.setAng(pair[0])
            self.setVel(pair[1])


    def setAng(self, ang):
        self.widgets['Shoot Angle'].set_val(ang)

    def setVel(self, vel):
        self.widgets['Shoot Speed'].set_val(vel)

    def updateAng(self, ang):
        if ang < -maxangle:
            ang = -maxangle
        elif ang > maxangle:
            ang = maxangle
        self.ang = ang
        self.go(run=False)

    def updateVel(self, vel):
        if vel < 0.01:
            print("Velocity must be >= 0.01")
            vel = 0.01
        self.vel = vel
        self.go(run=False)

    def run(self, tmax=None):
        self.model.compute('test', force=True)
        self.traj = self.model.trajectories['test']
        self.add_data_traj(self.traj)
        self.pts = self.points #Shouldn't have to do this.
        if self.calc_context is not None:
            # Update calc context
            self.calc_context()

    def get_forces(self, x, y):
        """
        For given x, y coord arguments, returns two dictionaries keyed
        by body number (1-N):
        net force magnitude, force vector
        """
        # Bombardier specific
        Fxs = []
        Fys = []
        Fs = []
        pars = self.model.query('pars')
        ixs = range(self.N)
        for i in ixs:
            m = pars['m%i'%i]
            bx = pars['bx%i'%i]
            by = pars['by%i'%i]
            p = pow(pp.distfun(x,y,bx,by),3)
            Fx = -m*(x-bx)/p
            Fy = -m*(y-by)/p
            Fxs.append(Fx)
            Fys.append(Fy)
            Fs.append(sqrt(Fx*Fx+Fy*Fy))
        return dict(zip(ixs, Fs)), dict(zip(ixs, zip(Fxs, Fys)))

    def set_planet_data(self, n, data):
        assert n in range(self.N)

        # default to old radius, unless updated (for masses)
        r = self.model.query('pars')['r%i'%n]
        d = self.density[n]
        pardict = {}
        for key, val in data.items():
            if key == 'r':
                pardict['r%i'%n] = val
                r = val
                self.radii[n] = r
            elif key == 'x':
                pardict['bx%i'%n] = val
                p = self.pos[n]
                self.pos[n] = (val, p.y)
            elif key == 'y':
                pardict['by%i'%n] = val
                p = self.pos[n]
                self.pos[n] = (p.x, val)
            elif key == 'd':
                d = val
                self.density[n] = d
            else:
                raise ValueError("Invalid parameter key: %s"%key)
            pardict['m%i'%n] = G*d*np.pi*r*r

        self.model.set(pars=pardict)
        self.body_pars.update(pardict)

        self.trajline = None
        self.startpt = None
        self.endpt = None
        self.go(run=False)
        self.graphics_refresh()

    def model_namer(self):
        name = 'sim_N%i'%self.N+'_fig%i'%self.fignum
        return name


# also make a circular version (using radial event)
class target4D_line(qt_feature_leaf):
    """
    Parameters expected:
    --------------------

    pt1, pt2 = Point2D specifying start and end of line in physical space
    speed, bearing = Interval objects for speed and bearing (may be singletons)
      N.B. bearing is measured relative to the perpendicular to the angle
           of the line, for convenience. I.e., 0 is fully transverse in direction
           of event detection direction, so typical intervals are [-45, 45]
    loc_event_name = name of zero-crossing event detecting goal location

    Assumptions:
    ------------
        System contains a uni-directional non-terminal event for the physical location

    """

    def evaluate(self, target):
        # target should be a model interface object
        ptsFS = target.test_traj.sample()
        pt1 = self.pars.pt1
        pt2 = self.pars.pt2
        speed_inter = self.pars.speed_inter
        bearing_inter = self.pars.bearing_inter
        ev_name = self.pars.loc_event_name

        # Compute metric for closeness along, and perpendicular to,
        # pt1-pt2 line.
        # Did zero crossing event occur with direction consistent with heading?
        event_ix_dict = ptsFS.labels.by_label['Event:'+ev_name]
        if event_ix_dict is not None:
            if len(event_ix_dict) != 1:
                raise ValueError
            else:
                ix = list(event_ix_dict.keys())[0]
                check_time = ptsFS['t'][ix]
                event_pt = ptsFS[ix]
        else:
            conditions_met = False
            # how close did the trajectory get, perpendicularly?
            print("Failed")
            return False

        # check that perp distance of event_pt to pt1-pt2 line is < epsilon
        #   tolerance of event
        # model = target.get('events', event_pt, check_time)
        ev = target.query('events')[ev_name]
        # initial value for conditions_met ...
        conditions_met = pp.distance_to_line(event_pt[['x','y']], (pt1, pt2)) < ev.eventtol # ?
        # get distance of event_pt (x,y) coords to pt1 and pt2
        # if either is > |pt1-pt2| then outside of sub-domain
        line_seg_len = np.linalg.norm(pt2 - pt1)
        dist1 = np.linalg.norm(event_pt[['x','y']] - pt1)
        dist2 = np.linalg.norm(event_pt[['x','y']] - pt2)
        conditions_met = conditions_met and ((dist1 < line_seg_len) and (dist2 < line_seg_len))
        # test if event_pt speed in speed_inter and bearing in bearing_inter
        conditions_met = conditions_met and event_pt['speed'] \
            in self.pars.speed_inter
        # define the sense of the line's angle based on event direction:
        # +1 perpendicular is angle 0
        line_vec = pt2 - pt1
        n = pp.get_orthonormal(line_vec)
        vel_vec = pp.Point2D({'x': event_pt['vx'], 'y': event_pt['vy']})
        speed = event_pt['speed']
        vel_vec = vel_vec/speed # normalized to unit length
        n_dot_v = np.dot(n, vel_vec)
        if n_dot_v < 0:
            # sense of n is backwards
            n = -n
        rel_bearing = pp.get_bearing(n, vel_vec)
        conditions_met = conditions_met and rel_bearing \
            in self.pars.bearing_inter
        # compute metric for closeness of velocity magnitude and direction
        return conditions_met



# ===========================================



def combine_planets(body_data, body1, body2):
    """
    Return effective combined single body from two bodies given
    as integers
    """
    r1 = body_data[body1]['radius']
    r2 = body_data[body2]['radius']
    m1 = np.pi*r1*r1
    m2 = np.pi*r2*r2
    d1 = body_data[body1]['density']
    d2 = body_data[body2]['density']
    pos1 = pp.Point2D(body_data[body1]['position'][0], body_data[body1]['position'][1])
    pos2 = pp.Point2D(body_data[body2]['position'][0], body_data[body2]['position'][1])
    dp = pos2 - pos1
    new_pos = pos2 - m1/(m1+m2)*dp
    # weighted mean of masses
    new_density = (m1*d1 + m2*d2)/(m1+m2)
    new_radius = float(sqrt((m1+m2)/new_density/np.pi))
    # mass is redundant
    return {'radius': new_radius, 'density': new_density,
            'mass': m1+m2, 'position': [new_pos[0], new_pos[1]]}


def make_forcelines(sim, x, y, i):
    #flines = []
    #for n in range(sim.N):
    #    i = n-1
    #Fsi = Fs[i]
    #Fvi_x, Fvi_y = Fvecs[i]
    bx, by = sim.pos[i]
    line_i = gx.line_GUI(sim, x, y, bx, by)
    return line_i

def relative_angle(line1, line2):
    # to help project line 1 onto line 2
    return line1.ang - line2.ang

def project(sim, x, y, i, pline):
    fline = make_forcelines(sim, x, y, i)
    Fs, Fvecs = sim.get_forces(x,y)
    F = Fs[i]
    ra = relative_angle(fline, pline)
    print(fline.ang_deg - pline.ang_deg)
    proj_F = F*cos(ra)
    bx, by = sim.pos[i]
    proj_d = pp.distfun(x,y,bx,by)*cos(ra)
    return proj_F, proj_d

def net_force_vs_line(sim, mesh_n=100):
    """
    Computes net forces of all bodies projected onto the line
    and to its normal direction.
    The line is assumed to be currently selected in the sim

    Mesh size given by optional argument mesh_n (default 100)
    """
    line = sim.selected_object
    if line is None:
        raise ValueError("No line selected")
    x1 = line.x1
    y1 = line.y1
    dx = line.dx
    dy = line.dy
    forces = np.zeros(mesh_n)
    forces_normal = np.zeros(mesh_n)
    for i, d in enumerate(np.linspace(0,1,mesh_n)):
        x = x1 + d*dx
        y = y1 + d*dy
        Fs, Fvecs = sim.get_forces(x,y)
        net_Fvec = np.sum(Fvecs.values(), 0)
        net_F = np.linalg.norm(net_Fvec)
        ang_Fvec = atan2(net_Fvec[1], net_Fvec[0])
        rel_ang = line.ang - ang_Fvec
        forces[i] = net_F*cos(rel_ang)
        forces_normal[i] = net_F*sin(rel_ang)
    return forces, forces_normal


def net_force_along_pts(sim, pts, bodies=None):
    """
    Raw force magnitudes along a pointset (not projected)
    for the selected bodies (by #), as well as total net force
    due to *all* bodies present

    If bodies is None (default) then all bodies with positive mass
    will be assumed.
    """
    num_pts = len(pts)
    if bodies is None:
        bodies = [i for i, d in enumerate(sim.density) if d > 0]
    Fs_by_body = {}
    net_Fs = np.zeros(num_pts)

    for n in bodies:
        Fs_by_body[n] = np.zeros(num_pts)

    for i, (x,y) in enumerate(pts[['x','y']].coordarray.T):
        Fs, Fvecs = sim.get_forces(x,y)
        for n in bodies:
            Fs_by_body[n][i] = Fs[n]
        try:
            net_Fs[i] = np.linalg.norm(np.sum(list(Fvecs.values()), 0))
        except TypeError:
            print(Fvecs.values())
            print(np.sum(list(Fvecs.values()), 0))

    return net_Fs, Fs_by_body


def plot_forces_along_line(sim, forces, fignum=2):
    line = sim.selected_object
    if line is None:
            raise ValueError("No line selected")
    plt.figure(fignum)
    mesh_n = len(forces)
    xaxis_range = linspace(0, line.length, mesh_n)
    plt.plot(xaxis_range, forces)
    plt.xlabel('distance along line')
    plt.ylabel('net force')
    plt.title('Forces vs. line %s' % line.name)
    plt.hlines(0, 0, line.length)
    axis('tight')


def pericenter_vs_n(sim, n, ecc, pt=None):
    # closest passing point relative to body n (index n-1)
    # at given point (or else IC assumed)
    if pt is None:
        pt = sim.model.query('ics')
    p0 = np.array((pt['x'], pt['y'])) # e.g. (0.25, 0.3)
    r0 = sim.pos[n]-p0
    v0  = np.array((pt['vx'], pt['vy']))
    v = np.linalg.norm(v0)
    r = np.linalg.norm(r0)
    m = sim.masses[n]
    mu = G*m
    # abs takes into account sign convention for hyperbolae
    # vs. ellipses
    eps = abs(v*v/2 - mu/r)
    a = mu / (2*eps)
    return abs((1-ecc)*a)

def apicenter_vs_n(sim, n, ecc, pt=None):
    # furthest point relative to body n (index n-1)
    # (assuming elliptical orbit, otherwise infinity is returned)
    # at given point (or else IC assumed)
    if pt is None:
        pt = sim.model.query('ics')
    if ecc >= 1:
        return np.Inf
    p0 = np.array((pt['x'], pt['y'])) # e.g. (0.25, 0.3)
    r0 = sim.pos[n]-p0
    v0  = np.array((pt['vx'], pt['vy']))
    v = np.linalg.norm(v0)
    r = np.linalg.norm(r0)
    m = sim.masses[n]
    mu = G*m
    # abs takes into account sign convention for hyperbolae
    # vs. ellipses
    eps = abs(v*v/2 - mu/r)
    a = mu / (2*eps)
    return abs((1+ecc)*a)

def eccentricity_vs_n(sim, n, pt=None):
    # relative to body n (index n-1)
    # at given point (or else IC assumed)
    if pt is None:
        pt = sim.model.query('ics')
    p0 = np.array((pt['x'], pt['y'])) # e.g. (0.25, 0.3)
    v0  = np.array((pt['vx'], pt['vy']))
    r0 = sim.pos[n]-p0
    r_cross_v0 = float(np.cross(r0, v0))
    m = sim.masses[n]
    mu = G*m
    v = np.linalg.norm(v0)
    r = np.linalg.norm(r0)
    e = sqrt(1+(v*v/(mu*mu) - 2/(r*mu))*(r_cross_v0)**2)
    return e


class body_context(fovea.calc_context):
    """
    Calculation context for a single planetary body with small satellite
    of negligible mass. n specifies which body.
    """
    def local_init(self, n, pt=None):
        self._refresh_init_args = (n, pt)
        if pt is None:
            pt = self.sim.model.query('ics')
        w = self.workspace
        w.pt = pt
        w.p0 = np.array((pt['x'], pt['y'])) # e.g. (0.25, 0.3)
        w.v0  = np.array((pt['vx'], pt['vy']))
        w.r0 = self.sim.pos[n]-w.p0
        # np.double might be easier to make work with Symbolic
        w.r_cross_v0 = np.double(np.cross(w.r0, w.v0))
        w.m = self.sim.masses[n]
        # ISSUE: G is a global!
        w.mu = G*w.m
        w.v = np.linalg.norm(w.v0)
        w.r = np.linalg.norm(w.r0)

class calc_context_forces(fovea.calc_context):
    def local_init(self):
        self.workspace.net_Fs, self.workspace.Fs_by_body = \
            net_force_along_pts(self.sim, self.sim.pts)

##def contextualize(context):
##    def decorator(target):
##        def wrapper():
##            return target(context)
##        return wrapper
##    return decorator
##
##@contextualize(my_context(game4, 1))
##def eccentricity(con):
##    return 1  # test
##
##def attach(con, fn):
##    def wrapped_fn():
##        return fn(con)
##    con.__dict__[fn.__name__] = wrapped_fn
##attach(con4, eccentricity)

@prep('ecc')
def eccentricity(con):
    return sqrt(1+(con.workspace.v*con.workspace.v/(con.workspace.mu*con.workspace.mu) \
                   - 2/(con.workspace.r*con.workspace.mu))*(con.workspace.r_cross_v0)**2)

@prep('eps')
def total_energy(con):
    # abs takes into account sign convention for hyperbolae
    # vs. ellipses
    return abs(con.workspace.v*con.workspace.v/2 - con.workspace.mu/con.workspace.r)

@prep('a')
def semimajor(con):
    return con.workspace.mu / (2*con.workspace.eps)

@prep('api')
def apicenter(con):
    return abs((1+con.workspace.ecc)*con.workspace.a)

@prep('peri')
def pericenter(con):
    return abs((1-con.workspace.ecc)*con.workspace.a)
