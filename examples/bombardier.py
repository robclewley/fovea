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
import fovea.domain2D as dom
from fovea import common

gentype = 'dopri'

if gentype == 'vode':
    targetlang = 'python'
else:
    targetlang = 'c'

# let rocket mass be negligible
# and G chosen to absorb m_rocket

# grav constant
G = 35

# Generic scale separation threshold
scale_thresh = 3.0

# graphics / control params
xdomain_halfwidth = .7
maxangle = 80 # degrees
da_dict = dict(zip( ('h','H','j','J','n','m'),
                    (-1, -10, 1, 10, -0.1, 0.1)))
dv_dict = dict(zip( ('d','D','f','F','c','v'),
                    (-1, -10, 1, 10, -0.1, 0.1)))

print "Change angle keys:"
print da_dict

print "Change velocity keys:"
print dv_dict

# other keys used in GUIrocket:
# l = make a line of interest (click-drag-release)
# SPACE = measure forces at clicked mouse point
# s = snap clicked mouse point to closest point on trajectory
# . (o) = grow a 2D domain
# g = GO! (run simulation)

dom_key = '.'
change_mouse_state_keys = ['l', 's', ' '] + [dom_key]

distfun = lambda x,y,x1,y1: sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))

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

class GUI_domain_handler(object):
    def __init__(self, gui, verbose=False):
        self.func = None
        self.domain = None  # polygon object
        self.gui = gui  # hook back
        self.polygon_domain_obj = None
        # grow states per polygon:
        #  0 = not begun, no domain object created
        #  1 = center placed, no initial radius
        #  2 = initial polygon defined, domain object created
        #  3 = final polygon defined
        #  4 = 3+1
        #  5 = 3+2
        #  (6 = new polygon defined) --> 3 after merge
        self.gui_grow_state = 0
        self.center_pt = None
        self.p1_pt = None
        self.verbose = verbose

    def assign_criterion_func(self, func):
        """
        func( (x,y) ) --> real scalar

        func is a zero-crossing function whose null-cline
        defines the boundary of a domain, or a sequence of
        such functions
        """
        self.func = func

    def event(self, evcode):
        """
        Mini finite state machine
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
        if self.verbose:
            print("In make_dom_c")
        # release mouse event control
        self.gui.fig.canvas.mpl_disconnect(self.gui.mouse_cid)
        # update state
        self.event('mouse')
        if self.gui_grow_state != 1:
            if self.verbose:
                print "make_dom_c failed"
            self.gui.mouse_wait_state_owner = None
            return
        # assign to c
        self.center_pt = pp.Point2D(ev.xdata, ev.ydata)
        # display c
        self.gui.selected_object_temphandle = self.gui.ax.plot(ev.xdata, ev.ydata, 'go')[0]
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
                print "make_dom_p1 failed"
            return
        # assign to p1
        self.p1_pt = pp.Point2D(ev.xdata, ev.ydata)
        # delete display of c
        self.gui.selected_object_temphandle.remove()
        # create then iterate polygon_domain object
        self.polygon_domain_obj = dom.polygon_domain(self.center_pt,
                                                     self.p1_pt,
                                                     self.func,
                                                     nsides=40,
                                                     edge_len_rtol=3)
        if self.verbose:
            print("Growing domain")
        self.polygon_domain_obj.grow()
        self.domain = self.polygon_domain_obj.polygon
        self.gui_grow_state = 3
        self.show_domain()
        if self.verbose:
            print("Domain complete")

    def show_domain(self):
        xs, ys = self.polygon_domain_obj.polygon.exterior.xy
##        if self.gui.selected_object_temphandle is not None:
            #print "show_domain ignoring: ", self.gui.selected_object_temphandle
##            try:
##                self.gui.selected_object_temphandle.remove()
##            except ValueError:
##                # sequence
##                for th in self.gui.selected_object_temphandle:
##                    th.remove()
        self.gui.selected_object_temphandle = self.gui.ax.plot(xs, ys, 'y-', lw=2, zorder=2)[0]
        self.gui.fig.canvas.draw()

    def unshow_domain(self):
        if self.gui.selected_object_temphandle is not None:
            try:
                self.gui.selected_object_temphandle.remove()
            except ValueError:
                # sequence
                for th in self.gui.selected_object_temphandle:
                    th.remove()
        self.gui.fig.canvas.draw()



class GUIrocket(object):
    def __init__(self, radii, density, pos, title, axisbgcol='black'):
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
        self.current_domain_handler = GUI_domain_handler(self)

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

        # Axes background colour
        self.axisbgcol = axisbgcol

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
        self.N = len(radii)
        self.gen_versioner_rocket = common.gen_versioner(os.path.abspath('.'),
                                                         self.name,
                                                         'simgen_N%i'%self.N,
                                                         gentype, 1)

        # Make this more generic for ABC
        self.setup_pars(radii, density, [pp.Point2D(x,y,labels={'body': i+1}) for i, (x,y) in enumerate(pos)])

        # --- END OF BOMBARDIER SPECIFICS

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

    def declare_in_context(self, con_obj):
        # context_changed flag set when new objects created and unset when Generator is
        # created with the new context code included
        self.context_changed = True
        self.context_objects.append(con_obj)

    def __str__(self):
        return self.name

    def setup_pars(self, radii, density, pos):
        # Should generalize to non-bombardier application
        N = self.N
        assert len(radii) == len(density) == len(pos) == N
        self.radii = radii
        self.density = density
        self.pos = pos # planet positions
        self.masses = [density[i]*np.pi*r*r for (i,r) in enumerate(radii)]
        rdict = dict([('r%i' %i, self.radii[i-1]) for i in range(1,N+1)])
        mdict = dict([('m%i' %i, self.masses[i-1]) for i in range(1,N+1)])
        posxdict = dict([('bx%i' %i, pos[i-1][0]) for i in range(1,N+1)])
        posydict = dict([('by%i' %i, pos[i-1][1]) for i in range(1,N+1)])
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
                self.model = self.gen_versioner_rocket.load_gen('sim_N%i'%self.N+'_fig%i'%self.fignum)
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
        for i in range(1,self.N+1):
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
            extra_auxvars.keys()
        DSargs.pars = pardict
        DSargs.pars.update(extra_pars)
        DSargs.fnspecs = auxfndict
        DSargs.fnspecs.update(extra_fnspecs)
        DSargs.algparams = {'init_step':0.001,
                            'max_step': 0.01,
                            'max_pts': 20000,
                            'maxevtpts': 2,
                            'refine': 5}

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
        for i in range(1,self.N+1):
            bev = Events.makeZeroCrossEvent('d(x,y,bx%i,by%i)-r%i' % (i,i,i),
                                            -1,
                                        {'name': 'b%iev' %i,
                                         'eventtol': 1e-5,
                                         'precise': True,
                                         'term': True},
                                        varnames=['x','y'],
                                        parnames=pardict.keys(),
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
        self.model = self.gen_versioner_rocket.make(DSargs)


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
        x = self.radii[-1]*cos(rad)
        y = -self.radii[-1]*sin(rad)
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
            print "Velocity must be >= 0.01"
            vel = 0.01
        self.vel = vel
        self.go(run=False)

    def key_on(self, ev):
        self._key = k = ev.key  # keep record of last keypress
        # TEMP
        print "Pressed", k
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
                self.ax.text(px-0.016,min(0.96, max(0.01,py-0.008)), str(i+1))

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
        self.model.compute('test')
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
        ixs = range(1,self.N+1)
        for i in ixs:
            m = pars['m%i'%i]
            bx = pars['bx%i'%i]
            by = pars['by%i'%i]
            p = pow(distfun(x,y,bx,by),3)
            Fx = -m*(x-bx)/p
            Fy = -m*(y-by)/p
            Fxs.append(Fx)
            Fys.append(Fy)
            Fs.append(sqrt(Fx*Fx+Fy*Fy))
        return dict(zip(ixs, Fs)), dict(zip(ixs, zip(Fxs, Fys)))

    def set_planet_data(self, n, data):
        assert n in range(1,self.N+1)

        # default to old radius, unless updated (for masses)
        r = self.model.query('pars')['r%i'%n]
        d = self.density[n-1]
        pardict = {}
        for key, val in data.items():
            if key == 'r':
                pardict['r%i'%n] = val
                r = val
                self.radii[n-1] = r
            elif key == 'x':
                pardict['bx%i'%n] = val
                p = self.pos[n-1]
                self.pos[n-1] = (val, p.y)
            elif key == 'y':
                pardict['by%i'%n] = val
                p = self.pos[n-1]
                self.pos[n-1] = (p.x, val)
            elif key == 'd':
                d = val
                self.density[n-1] = d
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
            data = findpt_nophase(self, ev.xdata, ev.ydata, eps=0.1)
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
        print("            = (%i, %.3f, (%.3f, %.3f))" % (data[0], data[1], x_snap, y_snap))
        self.fig.canvas.mpl_disconnect(self.mouse_cid)
        self.mouse_wait_state_owner = None

    def onselect_line(self, eclick, erelease):
        if eclick.button == 1:
            # left (primary)
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata
            self.selected_object = line_GUI(self, x1, y1, x2, y2)
            print("Created line as new selected object, now give it a name")
            print("  by writing this object's selected_object.name attribute")
            self.RS_line.set_active(False)
            self.mouse_wait_state_owner = None

    def onselect_box(self, eclick, erelease):
        self.mouse_wait_state_owner = None



def findpt_nophase(sim, x, y, eps=0.1, N_phases=20):
    try_phases = linspace(0, 1, N_phases)
    data = None
    min_d = 1000
    best_data = None
    for ph in try_phases:
        #print ph
        try:
            data = findpt(sim, x, y, ph, tol=1.0/N_phases, eps=eps)
        except ValueError:
            continue
        else:
            # found something
            if data[1] < min_d:
                best_data = data
                min_d = data[1]
    if best_data is None:
        raise ValueError("No point found")
    else:
        return best_data

def findpt(sim, x, y, phase, tol=0.2, eps=0.1):
    """
    Find closest point in pointset (not traj)
    by (x,y) and approximate phase in (0,1) fraction
    of whole trajectory.

    Optional tolerance fraction of #pts around phase
    guess to seek in (default 20%=0.2).

    Optional epsilon tolerance is the minimum distance
    from the desired point that target point must be.
    """
    pts = sim.pts
    n = len(pts)
    ixtol = int(tol*n)
    tguess = phase * pts['t'][-1]
    ixguess = pts.find(tguess, 0)
    ixlo = max(0, ixguess - ixtol)
    ixhi = min(n, ixguess + ixtol)
    new_pts = pts[ixlo:ixhi]
    ix1 = find(new_pts['x'], x)
    ix2 = find(new_pts['y'], y)
    d1 = distfun(new_pts[ix1]['x'], new_pts[ix1]['y'], x, y)
    d2 = distfun(new_pts[ix2]['x'], new_pts[ix2]['y'], x, y)
    if d1 < eps and d2 > eps:
        return ix1+ixlo, d1, new_pts[ix1]
    if d2 < eps and d1 > eps:
        return ix2+ixlo, d2, new_pts[ix2]
    if d1 < eps and d2 < eps:
        if d1 < d2:
            return ix1+ixlo, d1, new_pts[ix1]
        else:
            return ix2+ixlo, d2, new_pts[ix2]
    else:
        #print("findpt data was", (ix1, ix2), (d1, d2), (new_pts[ix1]['x'], new_pts[ix1]['y']),
        #      (new_pts[ix2]['x'], new_pts[ix2]['y']))
        #print '\n'
        raise ValueError("No match found to within eps=%.3f"%eps)


class line_GUI(object):
    def __init__(self, sim, x1,y1,x2,y2):
        if x1 > x2:
            # ensure correct ordering for angles
            xt = x1
            yt = y1
            x1 = x2
            y1 = y2
            x2 = xt
            y2 = yt
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.dy = y2-y1
        self.dx = x2-x1
        self.length = np.linalg.norm((self.dx, self.dy))
        # angle relative to horizontal, in radians
        self.ang = atan2(self.dy,self.dx)
        self.ang_deg = 180*self.ang/pi
        # hook back to sim object
        self.sim = sim
        # declare self to sim
        sim.declare_in_context(self)
        # move self to the currently selected object in sim
        sim.selected_object = self
        self.extra_fnspecs = {}
        self.extra_pars = {}
        self.extra_auxvars = {}
        self.extra_events = []
        print("Created line and moved to currently selected object")

        # actual MPL line object handle
        self.l = None
        self.name = '<untitled>'
        self.show()

    def __repr__(self):
        if self.extra_events == []:
            ev_str = '(no event)'
        else:
            ev_str = '(with event)'
        return "line_GUI(%.3f, %.3f, %.3f, %.3f) - '%s' %s" %(self.x1, self.y1, \
                                          self.x2, self.y2, self.name, ev_str)

    def show(self):
        if self.l is None:
            self.l = self.sim.ax.plot([self.x1,self.x2],
                                       [self.y1,self.y2],
                                   'y-')[0]
        else:
            self.l.set_visible(True)
        plt.draw()

    def unshow(self):
        self.l.set_visible(False)
        plt.draw()

    def remove(self):
        self.l.remove()
        plt.draw()

    def distance_to_pos(self, dist):
        """
        Calculate absolute (x,y) position of distance dist from (x1,y1) along line
        """
        return self.fraction_to_pos(self, dist/self.length)

    def fraction_to_pos(self, fraction):
        """
        Calculate absolute (x,y) position of fractional distance (0-1) from (x1,y1) along line
        """
        return np.array([self.x1+fraction*self.dx, self.y1+fraction*self.dy])

    def make_event_def(self, uniquename, dircode=0):
        self.name = uniquename
        res = pp.make_distance_to_line_auxfn('exit_line_'+uniquename, 'exit_fn_'+uniquename,
                                          ['x', 'y'], True)

        parname_base = 'exit_line_%s_' %uniquename
        self.extra_pars[parname_base+'p_x'] = self.x1
        self.extra_pars[parname_base+'p_y'] = self.y1
        self.extra_pars[parname_base+'dp_x'] = self.dx
        self.extra_pars[parname_base+'dp_y'] = self.dy
        self.extra_fnspecs.update(res['auxfn'])
        self.extra_events = [Events.makeZeroCrossEvent(expr='exit_fn_%s(x,y)' %uniquename,
                                              dircode=dircode,
                                              argDict={'name': 'exit_ev_%s' %uniquename,
                                                       'eventtol': 1e-8,
                                                       'eventdelay': 1e-3,
                                                       'starttime': 0,
                                                       'precise': True,
                                                       'active': True,
                                                       'term': False},
                                              varnames=('x', 'y'),
                                              fnspecs=res['auxfn'],
                                              parnames=res['pars'],
                                              targetlang=targetlang
                                              )]

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


class tracker_plotter(object):
    def __init__(self):
        self.figs = {}
        self.sim = None
        self.calc_context = None

    def __call__(self, calc_context, fignum, xstr, ystr, style):
        self.sim = calc_context.sim
        self.calc_context = calc_context
        fig = plt.figure(fignum)
        new_track = args(xstr=xstr, ystr=ystr, style=style)
        if fignum in self.figs:
            self.figs[fignum].tracked.append(new_track)
        else:
            self.figs[fignum] = args(figure=fig, tracked=[new_track])
        self.sim.tracked_objects.append(self)

    def show(self):
        for fignum, figdata in self.figs.items():
            plt.figure(fignum)
            ax = plt.gca()
            #figdata.figure.clf()
            ax.cla()
            for tracked in figdata.tracked:
                plt.plot(getattr(self.calc_context.workspace, tracked.xstr),
                     getattr(self.calc_context.workspace, tracked.ystr),
                     tracked.style, label=tracked.ystr.replace('_','\_'))
            plt.legend()
            plt.title('%s measures vs %s'%(self.calc_context.sim.name, tracked.xstr))
        plt.show()

# singleton
track_plot = tracker_plotter()



# ===========================================



def combine_planets(sim, body1, body2):
    """
    Return effective combined single body from two bodies given
    as integers
    """
    m1 = sim.masses[body1-1]
    m2 = sim.masses[body2-1]
    d1 = sim.density[body1-1]
    d2 = sim.density[body2-1]
    dp = sim.pos[body2-1] - sim.pos[body1-1]
    new_pos = sim.pos[body2-1] - m1/(m1+m2)*dp
    # weighted mean of masses
    new_density = (m1*d1 + m2*d2)/(m1+m2)
    new_radius = sqrt((m1+m2)/new_density/np.pi)
    return {'r': new_radius, 'd': new_density, 'm': m1+m2, 'p': new_pos}


def dist_vectorized(p1, p2vec, coordnames=None):
    """
    Vector of L2 distances between one point and one sequence of points,
    assumed to have same dimension unless optional tuple of coord names
    provided.
    """
    if coordnames is None:
        return np.array([np.linalg.norm(np.asarray(p1)-p) for p in p2vec])
    else:
        return np.array([np.linalg.norm(np.asarray(p1)-p) for p in p2vec[coordnames]])

def dist(p1, p2):
    """
    L2 distance between two points assumed to have same dimension
    """
    return np.linalg.norm(np.asarray(p2)-np.asarray(p1))


def make_forcelines(sim, x, y, i):
    #flines = []
    #for n in range(1,sim.N+1):
    #    i = n-1
    #Fsi = Fs[i]
    #Fvi_x, Fvi_y = Fvecs[i]
    bx, by = sim.pos[i-1]
    line_i = line_GUI(sim, x, y, bx, by)
    return line_i

def relative_angle(line1, line2):
    # to help project line 1 onto line 2
    return line1.ang - line2.ang

def project(sim, x, y, i, pline):
    fline = make_forcelines(sim, x, y, i)
    Fs, Fvecs = sim.get_forces(x,y)
    F = Fs[i]
    ra = relative_angle(fline, pline)
    print fline.ang_deg - pline.ang_deg
    proj_F = F*cos(ra)
    bx, by = sim.pos[i-1]
    proj_d = distfun(x,y,bx,by)*cos(ra)
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
        bodies = [i+1 for i, d in enumerate(sim.density) if d > 0]
    Fs_by_body = {}
    net_Fs = np.zeros(num_pts)

    for n in bodies:
        Fs_by_body[n] = np.zeros(num_pts)

    for i, (x,y) in enumerate(pts[['x','y']].coordarray.T):
        Fs, Fvecs = sim.get_forces(x,y)
        for n in bodies:
            Fs_by_body[n][i] = Fs[n]
        net_Fs[i] = np.linalg.norm(np.sum(Fvecs.values(), 0))

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


# Add this to utils.py or common.py in PyDSTool
def arclength(pts, vars=['x','y']):
    """
    Return array of L-2 arclength progress along parameterized pointset
    in the chosen dimensions
    """
    xpts = pts[vars]
    x0 = xpts[0]
    arclength = np.zeros(len(pts))
    for i, x in enumerate(xpts[1:]):
        arclength[i+1] = np.linalg.norm(x - xpts[i]) + arclength[i]
    return arclength


def pericenter_vs_n(sim, n, ecc, pt=None):
    # closest passing point relative to body n (index n-1)
    # at given point (or else IC assumed)
    if pt is None:
        pt = sim.model.query('ics')
    p0 = np.array((pt['x'], pt['y'])) # e.g. (0.25, 0.3)
    r0 = sim.pos[n-1]-p0
    v0  = np.array((pt['vx'], pt['vy']))
    v = np.linalg.norm(v0)
    r = np.linalg.norm(r0)
    m = sim.masses[n-1]
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
    r0 = sim.pos[n-1]-p0
    v0  = np.array((pt['vx'], pt['vy']))
    v = np.linalg.norm(v0)
    r = np.linalg.norm(r0)
    m = sim.masses[n-1]
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
    r0 = sim.pos[n-1]-p0
    r_cross_v0 = float(np.cross(r0, v0))
    m = sim.masses[n-1]
    mu = G*m
    v = np.linalg.norm(v0)
    r = np.linalg.norm(r0)
    e = sqrt(1+(v*v/(mu*mu) - 2/(r*mu))*(r_cross_v0)**2)
    return e


class calc_context(object):
    """
    __init__ method for concrete sub-class should insert any core parameters
    that are needed into 'shared' attribute
    """
    def __init__(self, sim, *args, **kwargs):
        self.sim = sim
        self._update_order = []
        self.workspace = dst.args()
        self._refresh_init_args = []
        # one function to one workspace variable
        self._functions_to_workspace = {}
        try:
            self.local_init(*args, **kwargs)
        except:
            print("local_init could not complete at initialization")

    def __call__(self):
        """
        Refresh workspace after update in attached simulator.
        Returns workspace.
        """
        self.local_init(*self._refresh_init_args)
        for fn_name in self._update_order:
            f = getattr(self, fn_name)
            # discard result but keep side-effects on workspace update
            f()
        return self.workspace


    def local_init(self, *args, **kwargs):
        """
        Optionally override in concrete sub-class
        """
        pass


    def attach(self, fn_seq):
        """Expect each function to have been decorated using
        @prep(<attr_name>)
        """
        if callable(fn_seq):
            # make a singleton list, for simplicity
            fn_seq = [fn_seq]
        for fn in fn_seq:
            self._attach(fn)


    def _attach(self, fn):
        """
        Seems that functions need to be wrapped individually
        in their own closure to avoid weird sharing of wrapped_fn
        """
        def wrapped_fn():
            val = fn(self)
            self.workspace[fn.attr_name] = val
            #print("Set workspace for value of %s is"%fn.attr_name, val)
            return val
        self.__setattr__(fn.__name__, wrapped_fn)
        self._functions_to_workspace[fn.__name__] = (wrapped_fn, fn.attr_name)

        # default to adding new function to end of update order
        self._update_order.append(fn.__name__)

        try:
            val = getattr(self, fn.__name__)() #fn(self)
        except:
            print("Could not compute value at attachment time for function %s"%fn.__name__)
             # initialize with None now, to declare in the meantime
            self.workspace[fn.attr_name] = None


class general_context(calc_context):
    """
    General purpose context
    """
    pass


class body_context(calc_context):
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
        w.r0 = self.sim.pos[n-1]-w.p0
        # np.double might be easier to make work with Symbolic
        w.r_cross_v0 = np.double(np.cross(w.r0, w.v0))
        w.m = self.sim.masses[n-1]
        w.mu = G*w.m
        w.v = np.linalg.norm(w.v0)
        w.r = np.linalg.norm(w.r0)


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

def prep(attr_name):
    def decorator(fn):
        fn.attr_name = attr_name
        return fn
    return decorator

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


def make_measure(fn_name, fn_spec, **defs):
    all_defs = defs.copy()
    q = dst.QuantSpec('_dummy_', fn_spec, treatMultiRefs=False)
    import PyDSTool.parseUtils as pu
    mapping = pu.symbolMapClass()
    assumed_modules = []
    tokens = q.parser.tokenized
    for sym in q.freeSymbols:
        # Hack, for now: if first (therefore, assumed all)
        # occurrences of symbol are in quotes, then don't convert.
        # Better solution would be to make parser create "x" as a single
        # symbol, at least with a detect quote option
        first_ix = tokens.index(sym)
        if first_ix == 0 or (first_ix > 0 and tokens[first_ix-1] not in ['"', "'"]):
            if pu.isHierarchicalName(sym):
                parts = sym.split('.')
                if parts[0] == 'sim':
                    mapping[sym] = 'con.'+sym
                elif parts[0] == 'bombardier':
                    # special case as this factory function is defined in that
                    # module so that reference will fail at runtime: remove
                    # 'bombardier' prefix
                    rest_sym = '.'.join(parts[1:])
                    mapping[sym] = rest_sym
                    scope = globals()
                    # locals override
                    scope.update(locals())
                    if parts[1] in scope:
                        all_defs[parts[1]] = scope[parts[1]]
                    else:
                        raise ValueError("Cannot resolve scope of symbol '%s'"%sym)
                else:
                    # assume module reference
                    assumed_modules.append(parts[0])
                    # record here to ensure inclusion in dyn_dummy
                    mapping[sym] = 'self.'+sym
            else:
                mapping[sym] = 'con.workspace.'+sym
        elif first_ix > 0 and tokens[first_ix-1] in ['"', "'"]:
            # put this symbol in the mapping values to ensure not included
            # as an argument to the function
            mapping[sym] = sym
    q.mapNames(mapping)
    import types
    for module_name in assumed_modules:
        global_scope = globals()
        # test if module name in scope
        if module_name in global_scope:
            _mod = global_scope[module_name]
            if isinstance(_mod, types.ModuleType):
                all_defs[module_name] = _mod

    # dyn_dummy contains dummy mappings but declares symbols to leave
    # evaluating until runtime
    dyn_dummy = dict(zip(mapping.values(), ['']*len(mapping)))
    funq = expr2fun(q, ensure_args=['con'], ensure_dynamic=dyn_dummy,
                   for_funcspec=False, fn_name=fn_name,
                   **all_defs)

    # decorate output
    funq.attr_name = fn_name
    return funq