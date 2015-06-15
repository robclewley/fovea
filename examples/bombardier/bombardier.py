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
from fovea import common, prep, graphics


gentype = 'dopri'

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

class GUIrocket(object):
    def __init__(self, bodies, title, axisbgcol='black'):
        """
        bodies is a dict-like mapping of 1 or more:
           <ID int>: {'density': <float>,
                      'radius': <float>,
                      'position': [<float>, <float>]}
        """
        global next_fighandle

        # Sim setup
        self.model = None
        # context objects (lines of interest, domains, etc)
        self.context_objects = []
        # external tracked objects (measures, etc)
        self.tracked_objects = []
        #
        self.selected_object = None
        self.selected_object_temphandle = None
        #
        self.current_domain_handler = dom.GUI_domain_handler(self)

        # name of task controlling mouse click event handler
        self.mouse_wait_state_owner = None

        # last output from a UI action
        self.last_output = None

        # if defined, will be refreshed on each Go!
        self.calc_context = None

        # Graphics widgets to be set for application
        self.widgets = {}


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
        # Currently unused
        #self.vtext = None
        #self.atext = None

        # ---------

        # ---- Make these generic in a parent class, then
        # specifcaly configured to bombardier here

        # Axes background colour
        self.axisbgcol = axisbgcol

        # Move these to a _recreate method than can be reused for un-pickling
        self.fig = figure(next_fighandle, figsize=(14,9))
        self.fignum = next_fighandle
        plt.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.1,
                               wspace=0.2, hspace=0.23)
        self.ax = plt.axes([0.05, 0.12, 0.9, 0.85], axisbg=axisbgcol)
        self.ax.set_title(title)
        self.name = title
        evKeyOn = self.fig.canvas.mpl_connect('key_press_event', self.key_on)
        evKeyOff = self.fig.canvas.mpl_connect('key_release_event', self.key_off)

        AngSlide = plt.axes([0.1, 0.055, 0.65, 0.03])
        self.widgets['AngBar'] = Slider(AngSlide, 'Shoot Angle', -maxangle, maxangle,
                                            valinit=self.ang, color='b',
                                            dragging=False, valfmt='%2.3f')
        self.widgets['AngBar'].on_changed(self.updateAng)

        VelSlide = plt.axes([0.1, 0.02, 0.65, 0.03])
        self.widgets['VelBar'] = Slider(VelSlide, 'Shoot Speed', 0.01, 2,
                                            valinit=self.vel, color='b',
                                            dragging=False, valfmt='%1.4f')
        self.widgets['VelBar'].on_changed(self.updateVel)

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
        GoButton = Button(plt.axes([0.005, 0.1, 0.045, 0.03]), 'Go!')
        GoButton.on_clicked(self.go)
        self.widgets['Go'] = GoButton

        self.RS_line = RectangleSelector(self.ax, self.onselect_line, drawtype='line') #,
#                                         minspanx=0.005, minspany=0.005)
        self.RS_line.set_active(False)
#        self.RS_box = RectangleSelector(self.ax, self.onselect_box, drawtype='box',
#                                        minspanx=0.005, minspany=0.005)
#        self.RS_box.set_active(False)

        # context_changed flag set when new objects created using declare_in_context(),
        # and unset when Generator is created with the new context code included
        self.context_changed = False
        self.setup_gen()
        self.traj = None
        self.pts = None

        self.mouse_cid = None # event-connection ID
        self.go(run=False)
        # force call to graphics_refresh because run=False above
        self.graphics_refresh(cla=False)
        plt.show()

        # next_fighandle for whenever a new model is put in a new figure (new game instance)
        next_fighandle += 1


    def graphics_refresh(self, cla=True):
        if cla:
            self.ax.cla()
        self.plot_bodies()
        self.plot_traj()
        # plot additional stuff
        self.plot_context()
        self.ax.set_aspect('equal')
        self.ax.set_xlim(-xdomain_halfwidth,xdomain_halfwidth)
        self.ax.set_ylim(0,1)
        self.fig.canvas.draw()

    def plot_context(self):
        for con_obj in self.context_objects:
            con_obj.show()
        for track_obj in self.tracked_objects:
            # external to main GUI window
            track_obj.show()

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

    def declare_in_context(self, con_obj):
        # context_changed flag set when new objects created and unset when Generator is
        # created with the new context code included
        self.context_changed = True
        self.context_objects.append(con_obj)

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

    def setup_gen(self):
        if self.context_changed:
            self.context_changed = False
            self.make_gen(self.body_pars, 'sim_N%i'%self.N+'_fig%i'%self.fignum)
        else:
            try:
                self.model = self.gen_versioner.load_gen('sim_N%i'%self.N+'_fig%i'%self.fignum)
            except:
                self.make_gen(self.body_pars, 'sim_N%i'%self.N+'_fig%i'%self.fignum)
            else:
                self.model.set(pars=self.body_pars)


    def make_gen(self, pardict, name):
        # scrape GUI diagnostic object extras for generator
        extra_events = []
        extra_fnspecs = {}
        extra_pars = {}
        extra_auxvars = {}
        for gui_obj in self.context_objects:
            extra_events.append(gui_obj.extra_events)
            extra_fnspecs.update(gui_obj.extra_fnspecs)
            extra_pars.update(gui_obj.extra_pars)
            extra_auxvars.update(gui_obj.extra_auxvars)

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
        self.fig.canvas.draw()
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
        self.widgets['AngBar'].set_val(ang)

    def setVel(self, vel):
        self.widgets['VelBar'].set_val(vel)

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

    def key_on(self, ev):
        self._key = k = ev.key  # keep record of last keypress
        # TEMP
        print("Pressed", k)
        if self.mouse_wait_state_owner == 'domain' and \
           k in change_mouse_state_keys:
            # reset state of domain handler first
            self.current_domain_handler.event('clear')

        if k in da_dict:
            Da = da_dict[k]*self.da
            self.updateAng(self.ang+Da)
            self.widgets['AngBar'].set_val(self.ang)
        elif k in dv_dict:
            Dv = dv_dict[k]*self.dv
            self.updateVel(self.vel+Dv)
            self.widgets['VelBar'].set_val(self.vel)
        elif k == 'g':
            print("Go! Running simulation.")
            self.go()
        elif k == 'l':
            print("Make a line of interest")
            self.RS_line.set_active(True)
            self.mouse_wait_state_owner = 'line'
        elif k == ' ':
            print("Forces at clicked mouse point")
            self.mouse_cid = self.fig.canvas.mpl_connect('button_release_event', self.mouse_event_force)
            self.mouse_wait_state_owner = 'forces'
        elif k == 's':
            print("Snap clicked mouse point to closest point on trajectory")
            self.mouse_cid = self.fig.canvas.mpl_connect('button_release_event', self.mouse_event_snap)
            self.mouse_wait_state_owner = 'snap'
        elif k == dom_key:
            print("Click on domain seed point then initial radius point")
            # grow domain
            if self.current_domain_handler.func is None:
                print("Assign a domain criterion function first!")
                return
            else:
                # this call may have side-effects
                self.current_domain_handler.event('key')
                self.mouse_wait_state_owner = 'domain'

    def key_off(self, ev):
        self._key = None

    def plot_bodies(self):
        for i in range(self.N):
            px, py = self.pos[i]
            if self.radii[i] > 0:
                if self.density[i] == 0:
                    col = 'green'
                else:
                    col = 'grey'
                self.ax.add_artist(plt.Circle((px,py),self.radii[i],color=col))
                self.ax.plot(px,py,'k.')
                self.ax.text(px-0.016,min(0.96, max(0.01,py-0.008)), str(i))

    def plot_traj(self, pts=None, with_speeds=True):
        """
        with_speeds option makes a "heat map" like color code along trajectory that denotes speed.
        """
        if pts is None:
            if self.pts is not None:
                pts = self.pts
            else:
                # nothing to plot
                return
        if self.axisbgcol == 'black':
            col = 'w'
        else:
            col = 'k'
        firstpt = pts[0]
        lastpt = pts[-1]

        if self.startpt is None:
            self.startpt = self.ax.plot(firstpt['x'],firstpt['y'],'ys', markersize=15)[0]
        else:
            self.startpt.set_xdata(firstpt['x'])
            self.startpt.set_ydata(firstpt['y'])

        if self.trajline is not None:
            self.trajline.remove()
        if with_speeds:
            speeds = pts['speed']
            norm = mpl.colors.Normalize(vmin=0, vmax=self.maxspeed)
            cmap=plt.cm.jet #gist_heat
            RGBAs = cmap(norm(speeds))
            xs = pts['x'][1:-1]
            ys = pts['y'][1:-1]
            segments = [( (xs[i], ys[i]), (xs[i+1], ys[i+1]) ) for i in range(len(xs)-1)]
            linecollection = mpl.collections.LineCollection(segments, colors=RGBAs)
            self.trajline = self.ax.add_collection(linecollection)
        else:
            self.trajline = self.ax.plot(pts['x'][1:-1], pts['y'][1:-1], col+'.-')[0]

        if self.endpt is None:
            self.endpt = self.ax.plot(lastpt['x'], lastpt['y'], 'r*', markersize=17)[0]
        else:
            self.endpt.set_xdata(lastpt['x'])
            self.endpt.set_ydata(lastpt['y'])
        n = len(pts)
        ptq1 = pts[int(0.25*n)]
        ptq2 = pts[int(0.5*n)]
        ptq3 = pts[int(0.75*n)]
        if self.quartiles is None:
            self.quartiles = [self.ax.plot(ptq1['x'], ptq1['y'], col+'d', markersize=10)[0],
                              self.ax.plot(ptq2['x'], ptq2['y'], col+'d', markersize=10)[0],
                              self.ax.plot(ptq3['x'], ptq3['y'], col+'d', markersize=10)[0]]
        else:
            self.quartiles[0].set_xdata(ptq1['x'])
            self.quartiles[0].set_ydata(ptq1['y'])
            self.quartiles[1].set_xdata(ptq2['x'])
            self.quartiles[1].set_ydata(ptq2['y'])
            self.quartiles[2].set_xdata(ptq3['x'])
            self.quartiles[2].set_ydata(ptq3['y'])
        plt.draw()


    def run(self, tmax=None):
        self.model.compute('test', force=True)
        self.traj = self.model.trajectories['test']
        self.pts = self.traj.sample()
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
        self.ax.cla()
        self.ax.set_aspect('equal')
        self.ax.set_xlim(-xdomain_halfwidth,xdomain_halfwidth)
        self.ax.set_ylim(0,1)

        self.trajline = None
        self.startpt = None
        self.endpt = None
        self.go(run=False)
        #self.graphics_refresh()

    def mouse_event_force(self, ev):
        print("\n(%.4f, %.4f)" %(ev.xdata, ev.ydata))
        fs, fvecs = self.get_forces(ev.xdata, ev.ydata)
        print(fs)
        print("Last output = (force mag dict, force vector dict)")
        self.last_output = (fs, fvecs)
        self.selected_object = pp.Point2D(ev.xdata, ev.ydata)
        if self.selected_object_temphandle is not None:
            self.selected_object_temphandle.remove()
        self.selected_object_temphandle = self.ax.plot(ev.xdata, ev.ydata, 'go')[0]
        self.fig.canvas.draw()
        self.fig.canvas.mpl_disconnect(self.mouse_cid)
        self.mouse_wait_state_owner = None

    def mouse_event_snap(self, ev):
        if self.pts is None:
            print("No trajectory defined")
            return
        print("\nClick: (%.4f, %.4f)" %(ev.xdata, ev.ydata))
        # have to guess phase, use widest tolerance
        try:
            data = pp.find_pt_nophase_2D(self.pts, pp.Point2D(ev.xdata, ev.ydata),
                                         eps=0.1)
        except ValueError:
            print("No nearby point found. Try again")
            self.fig.canvas.mpl_disconnect(self.mouse_cid)
            return
        self.last_output = data
        x_snap = data[2]['x']
        y_snap = data[2]['y']
        self.selected_object = pp.Point2D(x_snap, y_snap)
        if self.selected_object_temphandle is not None:
            self.selected_object_temphandle.remove()
        self.selected_object_temphandle = self.ax.plot(x_snap, y_snap, 'go')[0]
        self.fig.canvas.draw()
        print("Last output = (index, distance, point)")
        print("            = (%i, %.3f, (%.3f, %.3f))" % (data[0], data[1],
                                                          x_snap, y_snap))
        self.fig.canvas.mpl_disconnect(self.mouse_cid)
        self.mouse_wait_state_owner = None

    def onselect_line(self, eclick, erelease):
        if eclick.button == 1:
            # left (primary)
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata
            self.selected_object = graphics.line_GUI(self, self.ax,
                                            pp.Point2D(x1, y1),
                                            pp.Point2D(x2, y2))
            print("Created line as new selected object, now give it a name")
            print("  by writing this object's selected_object.name attribute")
            self.RS_line.set_active(False)
            self.mouse_wait_state_owner = None

    def onselect_box(self, eclick, erelease):
        self.mouse_wait_state_owner = None





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
    line_i = graphics.line_GUI(sim, x, y, bx, by)
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
