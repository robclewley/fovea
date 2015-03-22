"""
Scenario 1

Game target = source

General situation (full game 1)
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



Single body (reduced game 4)
----------------------------

Analytically solvable case!

This virtual body is offset from center of space

This sub-problem has a greatly regularized fitness landscape

Feature requirements:
    Put half-way point opposite body
    Ensure pericenter + apicenter < length of domain
      and apicenter close to initial distance from body
      Metric: |distance to body - apicenter| (minimize)


"""
from __future__ import division
from PyDSTool import *
import PyDSTool.Toolbox.phaseplane as pp
from bombardier import *
import bombardier


class calc_context_forces(bombardier.calc_context):
    def local_init(self):
        self.workspace.net_Fs, self.workspace.Fs_by_body = \
            bombardier.net_force_along_pts(self.sim, self.sim.pts)


# -------------------------------------

# Scenario specification
# ! W1a Objects:
density1 = [0.5, 0.5, 0.5, 0.5, 0.0]
radii1 = [0.04, 0.04, 0.03, 0.03, 0.08]
pos1 = [(-0.2, 0.45), (0.2, 0.7), (0.3, 0.4), (-0.14,0.55),
        (0,0)]

game1 = GUIrocket(radii1, density1, pos1, "Demo-game1", axisbgcol='white')
# ! W1b Initial conditions
game1.set( (-79, 0.7) )

# ! W2a Constraints
# NONE

# ! W2b Goals / targets
# User interaction to draw line

#ltarget = game1.selected_object

ltarget = line_GUI(game1, 0.36, 0.74, 0.42, 0.8)
ltarget.make_event_def('target1', 1)
game1.setup_gen()

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

game1.go()
test_model = intModelInterface(game1.model)
print("Success? %s"%(str(target(test_model))))
1/0
# ! W3 Create exploratory tool & calculation context objectsc

##game1.calc_context = calc_context_forces(game1)
##con1 = game1.calc_context
##w1 = con1.workspace
##uniformity_force = make_measure('uniformity_force', 'math.sqrt(np.std(net_Fs))')
##con1.attach(uniformity_force)

# alternative (deprecated) method
ecc1 = eccentricity_vs_n(game1, 1)
peri1 = pericenter_vs_n(game1, 1, ecc1)

1/0

# ============  SCENARIO 2  ============

# hypothesize combination of bodies 1 & 4 and 2 & 3
# then the combo of those
# provided stay sufficiently far from them
reduced_data14 = combine_planets(game1, 1, 4)
reduced_data23 = combine_planets(game1, 2, 3)

density2 = [reduced_data14['d'], reduced_data23['d'], 0.0]
radii2 = [reduced_data14['r'], reduced_data23['r'], 0.08]
pos2 = [reduced_data14['p'], reduced_data23['p'], pos1[4]]
game2 = GUIrocket(radii2, density2, pos2, "Demo-reduced-game2", axisbgcol='black')
game2.set( (-70, 0.7) )

game2.go()

uniformity_force = make_measure('uniformity_force', 'math.sqrt(np.std(net_Fs))')

arc = make_measure('arclength', 'bombardier.arclength(sim.pts)')

contrib_1to2 = make_measure('contrib_1to2', 'Fs_by_body[1]/Fs_by_body[2]')

contrib_1toall = make_measure('contrib_1toall', 'Fs_by_body[1]/net_Fs')

dom_thresh = 0.6

def body1_dominant_at_point(pt):
    """
    Returns scalar relative to user threshold of %age dominant out of net force
    """
    global dom_thresh
    net_Fs = game2.get_forces(pt[0],pt[1])[0]
    return net_Fs[1]/sum(net_Fs.values()) - dom_thresh

def body4_dominant_at_point(pt):
    """
    Returns scalar relative to user threshold of %age dominant out of net force
    """
    global dom_thresh
    net_Fs = game1.get_forces(pt[0],pt[1])[0]
    return net_Fs[4]/sum(net_Fs.values()) - dom_thresh

game1.current_domain_handler.assign_criterion_func(body4_dominant_at_point)

# Domain growth testing
##game2.current_domain_handler.assign_criterion_func(body1_dominant_at_point)
##1/0
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

# scalar measure of relative uniformity of net force over entire orbit
# (surrogate for eccentricity in multi-body system)
#@prep('uniformity_force')
#def uniformity_force(con):
#    return np.std(con.workspace.net_Fs)

# Now, want to attach these measures to a context and select which to be hooked up
# to auto-updated plots as game2 is refreshed
game2.calc_context = calc_context_forces(game2)
con2 = game2.calc_context
w2 = con2.workspace

con2.attach((arc, contrib_1to2, contrib_1toall,
                           uniformity_force))
#con2._update_order = ['', '']

# skip track plot demo here

##track_plot(con2, 10, 'arclength', 'contrib_1to2', 'k--')
##track_plot(con2, 10, 'arclength', 'contrib_1toall', 'k:')
##track_plot.show()

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

#eccentricity = make_measure('ecc', 'sqrt(1+(v*v/(mu*mu) - 2/(r*mu))*(r_cross_v0)**2)')
#context.attach(eccentricity)

# This calc_context only applies when body 1 gravity is suficiently strong
con2_vs1 = bombardier.body_context(game2, 1)
con2_vs1.attach((eccentricity, total_energy, semimajor, apicenter, pericenter))
w = con2_vs1.workspace

# decide to monitor actual distance from body 1 along orbit
# to compare with analytical estimate
dist_to_1_vs_peri = make_measure('dist_to_1_vs_peri', 'bombardier.dist_vectorized(sim.pos[0], sim.pts[["x","y"]]) - workspace1.peri', workspace1=w)
con2.attach(dist_to_1_vs_peri)

track_plot(con2, 11, 'arclength', 'dist_to_1_vs_peri', 'k-')
track_plot.show()
print "Mismatch of pericenter prediction without reduction", abs(min(w2.dist_to_1_vs_peri))
game2.current_domain_handler.assign_criterion_func(body1_dominant_at_point)

"""
We find the comparison of analytic peri based on IC to actual peri is not good,
presumably because initial condition is not strongly in the confidence zone
for body 1's dominance. Need error estimates and correction using Game 4 (single combined body).
"""

# skip Game 3 for now, which has one more regular planet restored

### ============  GAME 3  ============
##
##ecc2 = eccentricity_vs_n(game2, 1)
##peri2 = pericenter_vs_n(game2, 1, ecc2)
##
##
##density3 = [reduced_data14['d'], density1[1], density1[2], 0.0]
##radii3 = [reduced_data14['r'], radii1[1], radii1[2], 0.08]
##pos3 = [reduced_data14['p'], pos1[1], pos1[2], pos1[4]]
##game3 = GUIrocket(radii3, density3, pos3, "Demo-reduced-game3", axisbgcol='black')
##game3.set( (-72, 0.7) )
##
##ecc3 = eccentricity_vs_n(game3, 1)
##peri3 = pericenter_vs_n(game3, 1, ecc3)


# ============  GAME 4  ============

reduced_data_all = combine_planets(game2, 1, 2)

density4 = [reduced_data_all['d'], 0.0]
radii4 = [reduced_data_all['r'], 0.08]
pos4 = [reduced_data_all['p'], pos1[4]]
game4 = GUIrocket(radii4, density4, pos4, "Demo-reduced-game4", axisbgcol='black')
game4.set( (-72, 0.7) )
game4.go()


game4.calc_context = bombardier.body_context(game4, 1)
con4 = game4.calc_context
con4.attach((eccentricity, total_energy, semimajor, apicenter, pericenter))
w = con4()  # refreshes calcs, returns workspace
ecc4 = w.ecc
peri4 = w.peri
api4 = w.api

game4.set( (-70, 0.7) )
game4.go()

1/0

"""
Game 2:
  Define body-1-dominant region
  Define line at the edge of that region

Game 4:
  Assign that line to game4's currently selected object
  Name it
  Create the event for it
  Compute end point of this regime

Game 2:
  Assign new IC from end point of Game 4 regime and predict whether pericenter
    estimate is now accurate
"""

l2 = game2.selected_object
l4 = line_GUI(game4, l2.x1, l2.y1, l2.x2, l2.y2)
l4.name = 'regime1_bd'
l4.make_event_def('reg1_bd', 1)
game4.setup_gen()
game4.go()
exit_pt = game4.traj.getEvents('exit_ev_reg1_bd')[0]
game2.set(exit_pt[['vx','vy']], exit_pt[['x','y']], by_vel=True)
game2.run()
game2.graphics_refresh(cla=False)
print "Mismatch of pericenter prediction with reduction", abs(min(w2.dist_to_1_vs_peri))