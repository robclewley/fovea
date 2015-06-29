"""
Scenario 1: Game 2

Two bodies (reduced game)
-------------------------

Hypothesize combination of bodies 1 & 4 and 2 & 3
provided trajectory stays sufficiently far from them

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

body_setup1 = setup['Full_model']

# hypothesize combination of bodies 1 & 4 and 2 & 3
# then the combo of those
# provided stay sufficiently far from them
reduced_data_14 = {1: combine_planets(body_setup1, 1, 4)}
reduced_data_23 = {2: combine_planets(body_setup1, 2, 3)}

body_setup2 = {}
body_setup2.update({0: body_setup1[0]})
body_setup2.update(reduced_data_14)
body_setup2.update(reduced_data_23)

with open('bodies_setup.yaml', 'a') as f:
    yaml.dump({'Reduced_model2': body_setup2}, f)

game2 = GUIrocket(body_setup2, "Scenario 1: Game 2", axisbgcol='black')
game2.set( (-70, 0.7) )

game2.go()

dom_thresh = 0.6

def body1_dominant_at_point(pt_array, fsign=None):
    """
    Returns scalar relative to user threshold of %age dominant out of net force
    """
    global dom_thresh
    net_Fs = game2.get_forces(pt_array[0],pt_array[1])[0]
    return net_Fs[1]/sum(list(net_Fs.values())) - dom_thresh


# Domain growth testing
game2.current_domain_handler.assign_criterion_func(body1_dominant_at_point)

##
##pd = game2.current_domain_handler.polygon_domain_obj
##n=0
##import sys
##def grow_step():
##    global n
##    n+=1
##    print(n)
##    sys.stdout.flush()
##    pd.grow_step(verbose=True)
##    game2.current_domain_handler.show_domain()
##    if any(pd.stop_growing):
##        print(pd.stop_growing)
##        sys.stdout.flush()

# -------------------------------------------------
# Deprecated ways of obtaining measures of interest

#@prep('arclength')
#def arc(con):
#    # vector of arclength along pts
#    return arclength(con.sim.pts)

#@prep('contrib_1to2')
#def contrib_1to2(con):
#    # vector measure of relative contribution by body 1 along orbit
#    return con.workspace.Fs_by_body[1]/con.workspace.Fs_by_body[2]

#@prep('contrib_1toall')
#def contrib_1toall(con):
#    return con.workspace.Fs_by_body[1]/con.workspace.net_Fs

# scalar measure of relative variability of net force over entire orbit
# (surrogate for eccentricity in multi-body system)
#@prep('variability_force')
#def variability_force(con):
#    return np.std(con.workspace.net_Fs)
# -------------------------------------------------

# Attach measures to a context and select which to be hooked up
# to auto-updated plots as game2 is refreshed
game2.calc_context = calc_context_forces(game2, 'con2')
con2 = game2.calc_context

variability_force = fovea.make_measure('variability_force',
                                 'math.sqrt(np.std(net_Fs))')

con2.declare('PyDSTool.Toolbox.phaseplane', 'pp')
arc = fovea.make_measure('arclength', 'pp.arclength(sim.pts)')

contrib_1to2 = fovea.make_measure('contrib_1to2',
                            'Fs_by_body[1]/Fs_by_body[2]')

contrib_1toall = fovea.make_measure('contrib_1toall',
                              'Fs_by_body[1]/net_Fs')

w2 = con2.workspace


con2.attach((arc, contrib_1to2, contrib_1toall,
                           variability_force))
# Don't need to update order
#con2._update_order = ['', '']

#print(con2.arclength())
#print(con2.contrib_1to2())

tracker(con2, 10, 'arclength', 'contrib_1to2', 'k--')
tracker(con2, 10, 'arclength', 'contrib_1toall', 'k:')
tracker.show()

"""
Task: Find zone in which Body 1 dominates

Define by contrib_1toall > 0.75

1. Select good initial position  in zone
2. Grow zone as a circle until fails
(3.) Grow similar circular zones and merge them -> tubes, etc.
4. Create bounding events

... Why? It establishes zone in which analytical estimates
are most valid. So, if we start ICs in that zone then we don't
need simulations to estimate outcomes at exit point

How might it help generate or use further reduced models?
   Generate: analytical estimates or a single-body model apply within domain
   Error bounds / estimates / correction due to single-body approximation?
How might it help gain insight into the original, full model?
   Enables intuitive partial trajectory plan through this zone,
      and how it depends on the control parameters a priori

Task: Find sensitivity of zone boundaries to system uncertainty
Q. How does the zone change as system params (masses, positions, etc.) vary?
   E.g. to explore their uncertainty to find safest trajectory area to use


"""

eccentricity = fovea.make_measure('ecc', 'sqrt(1+(v*v/(mu*mu) - 2/(r*mu))*(r_cross_v0)**2)')

# This calc_context only applies when body 1 gravity is sufficiently strong
con2_vs1 = bombardier.body_context(game2, 'con2_vs1', 1)
con2_vs1.declare('PyDSTool.Toolbox.phaseplane', 'pp')

# Reuse existing arc measure object
# arc = fovea.make_measure('arclength', 'pp.arclength(sim.pts)')

con2_vs1.attach((arc, eccentricity, total_energy, semimajor, apicenter,
                 pericenter))
w = con2_vs1.workspace

# decide to monitor actual distance from body 1 along orbit
# to compare with analytical estimate
dist_to_1_vs_peri = fovea.make_measure('dist_to_1_vs_peri', 'pp.dist_vectorized(sim.pos[0], sim.pts[["x","y"]]) - workspace1.peri', workspace1=w)
con2_vs1.attach(dist_to_1_vs_peri)

tracker(con2_vs1, 11, 'arclength', 'dist_to_1_vs_peri', 'k-')
tracker.show()
print("Mismatch of pericenter prediction without reduction: %.5f" % abs(min(w.dist_to_1_vs_peri)))
game2.current_domain_handler.assign_criterion_func(body1_dominant_at_point)

# === FAILS!
#saveObjects(game2, 'game2.sav', force=True)

"""
We find the comparison of analytic peri based on IC to actual peri is not good,
presumably because initial condition is not strongly in the confidence zone
for body 1's dominance. Need error estimates and correction using Game 4 (single combined body).
"""
