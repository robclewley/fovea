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

import yaml
with open('bodies_setup.yaml') as f:
    setup = yaml.load(f)

# -------------------------------------

# Scenario specification (codes refer to usage document)
# ! W1a Objects:
body_setup1 = setup['Full_model']

game1 = GUIrocket(body_setup1, "Scenario 1: Game 1", axisbgcol='white')
# ! W1b Initial conditions
game1.set( (-79, 0.7) )

# ! W2a Constraints
# NONE

# ! W2b Goals / targets
# User interaction to draw line

#ltarget = game1.selected_object

ltarget = gx.line_GUI(game1, pp.Point2D(0.36, 0.74),
                      pp.Point2D(0.42, 0.8), subplot='11')
ltarget.make_event_def('target1', 1)
game1.setup_gen(game1.model_namer)

# make event terminal
game1.model.setDSEventTerm('gen', 'exit_ev_target1', True)

target = target4D_line('test_line', pars=args(pt1=pp.Point2D((ltarget.x1, ltarget.y1)),
                                          pt2=pp.Point2D((ltarget.x2, ltarget.y2)),
                                          speed_inter=Interval('speed', float, (1,2)),
                                          bearing_inter=Interval('bearing', float, (-15,45)),
                                          loc_event_name='exit_ev_target1')
                       )

#l = array((target.pars.pt1, target.pars.pt2))
#game1.ax.plot(l.T[0], l.T[1], 'k-', lw=3)
#plt.draw()

# ! W3 Create exploratory tool & calculation context objects

game1.calc_context = calc_context_forces(game1, 'con1')
con1 = game1.calc_context
w1 = con1.workspace
variability_force = fovea.make_measure('variability_force', 'math.sqrt(np.std(net_Fs))')
con1.attach(variability_force)

game1.go()
test_model = intModelInterface(game1.model)
print("Success? %s"%(str(target(test_model))))

print("Variability of net force felt along trajectory = %.3f" % con1.variability_force())
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
    net_Fs = game1.get_forces(pt_array[0],pt_array[1])[0]
    return net_Fs[4]/sum(list(net_Fs.values())) - dom_thresh

game1.current_domain_handler.assign_criterion_func(body4_dominant_at_point)

"""
User now clicks near body 4 to center the initial domain, and again a little further
away to seed the initial radius.
"""
1/0

# Run this to show domain growth.
pd = game1.current_domain_handler.polygon_domain_obj
n=0
import sys
def grow_step():
    global n
    n+=1
    print(n)
    sys.stdout.flush()
    pd.grow_step(verbose=True)
    game1.current_domain_handler.show_domain()
    if any(pd.stop_growing):
        print(pd.stop_growing)
        sys.stdout.flush()

