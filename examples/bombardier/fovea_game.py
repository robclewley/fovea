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
import fovea
import fovea.graphics as gx
from fovea.graphics import tracker

import fovea.domain2D as dom
from fovea import common, prep

import yaml
with open('bodies_setup.yaml') as f:
    setup = yaml.load(f)

# Scenario specification (codes refer to usage document)
# ! W1a Objects:
body_setup1 = setup['Reduced_model4']

game1 = GUIrocket(body_setup1, "Scenario 1: Game 4", axisbgcol='white')
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

##game1.calc_context = calc_context_forces(game1, 'con1')
##con1 = game1.calc_context
##w1 = con1.workspace
##variability_force = fovea.make_measure('variability_force', 'math.sqrt(np.std(net_Fs))')
##con1.attach(variability_force)

con = body_context(game1, 'con', 1)
con.attach((eccentricity, total_energy, semimajor, apicenter,
                 pericenter))

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