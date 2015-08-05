"""
Scenario 1: Game 1

General situation (full game)
-------------------------------

General feature requirements:
    Convex trajectory (no points of inflection) i.e. curvature > 0
      Metric: curvature
    Eccentricity < 1 (elliptical or circular orbit)
      Metric: eccentricity around a body
    Source is on the ellipse or circle
      Event detection determines
    Trajectory must not hit boundaries
      Non-detection of boundary events
    Trajectory stays far from bodies


Test assumption:
    Calculate eccentricity using distance of #1 from source
      as a (near-circular) radius

Test optimization:
    Calculate bearing and speed that will achieve low
      eccentricty

"""
from __future__ import division
from PyDSTool import *
import PyDSTool.Toolbox.phaseplane as pp
from bombardier import *
import bombardier
import fovea
import fovea.graphics as gx
from fovea.graphics import tracker

import fovea.domain2D as dom
from fovea import common, prep, graphics

import yaml
with open('bodies_setup.yaml') as f:
    setup = yaml.load(f)

# -------------------------------------
class GUIrocket(gx.diagnosticGUI):
    def __init__(self, bodies, title, axisbgcol='black'):
        """
        bodies is a dict-like mapping of 1 or more:
           <ID int>: {'density': <float>,
                      'radius': <float>,
                      'position': [<float>, <float>]}
        """
        global next_fighandle

        plotter = gx.plotter2D()
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
        self.addFig('master',
                        title='Bombardier',
                        xlabel='x', ylabel='y',
                        domain=DOI)

        #Setup all layers

        self.addLayer('trajs')
        self.addLayer('bodies', kind='patch')
        self.addLayer('text', kind='text')

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

        fig_struct, fig = self.plotter._resolveFig('master')
        self.ax = fig_struct.arrange['11']['axes_obj']

        self.addWidget(Slider, callback=self.updateAng, axlims = (0.1, 0.055, 0.65, 0.03),
                      label='Shoot Angle', valmin= -maxangle, valmax= maxangle,
                      valinit= self.ang, color='b', dragging=False, valfmt='%2.3f')

        self.addWidget(Slider, callback=self.updateVel, axlims=(0.1, 0.02, 0.65, 0.03),
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

        self.addWidget(Button, callback=self.go, axlims=(0.005, 0.1, 0.045, 0.03), label='Go!')

        # context_changed flag set when new objects created using declare_in_context(),
        # and unset when Generator is created with the new context code included
        self.context_changed = False
        self.setup_gen(self.model_namer)

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
            self.addDataPoints(quarts, coorddict=coorddict)

        except TypeError:
            pass

        #Traj Pointset
        coorddict = {'x':
                     {'x':'x', 'y':'y','layer':'trajs','name':'data1', 'object':'collection'},
                      #{'x':'x', 'y':'y','layer':'trajs', 'object':'collection'},
                     'speed':
                     {'map_color_to':'x'}
                     }
        self.addDataPoints(self.points, coorddict=coorddict)

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
        self.addDataPoints(bodsPoints, coorddict=coorddict)

        pos = np.array(self.pos).transpose()
        for i in range(len(pos[0])):
            self.plotter.addText(pos[0][i], pos[1][i], i, style='k', layer='text')

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
        self.addDataTraj(self.traj)
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

    def onselect_box(self, eclick, erelease):
        self.mouse_wait_state_owner = None

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

    def model_namer(self):
        name = 'sim_N%i'%self.N+'_fig%i'%self.fignum
        return name


# Scenario specification (codes refer to usage document)
# ! W1a Objects:
body_setup1 = setup['Full_model']

game1 = GUIrocket(body_setup1, "Scenario 1: Game 1", axisbgcol='white')
# ! W1b Initial conditions
game1.set( (-79, 0.7) )

ltarget = gx.line_GUI(game1, pp.Point2D(0.36, 0.74),
                      pp.Point2D(0.42, 0.8), subplot = '11')

ltarget.make_event_def('target1', 1)
game1.setup_gen(game1.model_namer)

## make event terminal
game1.model.setDSEventTerm('gen', 'exit_ev_target1', True)

target = target4D_line('test_line', pars=args(pt1=pp.Point2D((ltarget.x1, ltarget.y1)),
                                          pt2=pp.Point2D((ltarget.x2, ltarget.y2)),
                                          speed_inter=Interval('speed', float, (1,2)),
                                          bearing_inter=Interval('bearing', float, (-15,45)),
                                          loc_event_name='exit_ev_target1'))

game1.calc_context = calc_context_forces(game1, 'con1')
con1 = game1.calc_context
w1 = con1.workspace
variability_force = fovea.make_measure('variability_force', 'math.sqrt(np.std(net_Fs))')
con1.attach(variability_force)

game1.go()
test_model = intModelInterface(game1.model)
print("Success? %s"%(str(target(test_model))))

#print("Variability of net force felt along trajectory = %.3f" % con1.variability_force())
print(" (smaller is better)")

# alternative (deprecated) method, shown for reference
ecc1 = eccentricity_vs_n(game1, 1)
peri1 = pericenter_vs_n(game1, 1, ecc1)
print("Eccentricity = %.3f" % ecc1)

dom_thresh = 0.6

def body4_dominant_at_point(pt_array, fsign=None):
    """
    Returns scalar relative to user threshold of %age dominant out of net force
    """
    global dom_thresh
    net_Fs = game1.user_func(pt_array[0],pt_array[1])[0]
    return net_Fs[4]/sum(list(net_Fs.values())) - dom_thresh

game1.assign_user_func(game1.get_forces)
game1.current_domain_handler.assign_criterion_func(body4_dominant_at_point)

fig_struct, figure = game1.plotter._resolveFig(None)

#class snap_point():
    #def __init__(self, x, y):
        #self.xdata = x
        #self.ydata = y

#sp = snap_point(-0.51, 0.51)
#game1.mouse_event_snap(sp)

halt = True