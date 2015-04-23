"""
Scenario 1: Game 4

Single body (reduced game)
----------------------------

Combination of both reduced bodies from Game 2

Analytically solvable case!

This virtual body is offset from center of space

This sub-problem has a greatly regularized fitness landscape

Feature requirements:
    Put half-way point opposite body
    Ensure pericenter + apicenter < length of domain
      and apicenter close to initial distance from body
      Metric: |distance to body - apicenter| (minimize)

"""
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

body_setup2 = setup['Reduced_model2']
reduced_data_all = {1: combine_planets(body_setup2, 1, 2)}

body_setup4 = {}
body_setup4.update({0: body_setup2[0]})
body_setup4.update(reduced_data_all)

with open('bodies_setup.yaml', 'a') as f:
    yaml.dump({'Reduced_model4': body_setup4}, f)

game4 = GUIrocket(body_setup4, "Scenario 1: Game 4", axisbgcol='black')
game4.set( (-72, 0.7) )
game4.go()

game4.calc_context = bombardier.body_context(game4, "con4", 1)
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
Interactive scenario continuation:

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

# LOAD game2! (not yet implemented)
#l2 = game2.selected_object

# Create l2 graphically, based on data from a concurrent session with game2
l4 = gx.line_GUI(game4, game4.ax, pp.Point2D(l2.x1, l2.y1),
                 pp.Point2D(l2.x2, l2.y2))
l4.name = 'regime1_bd'
l4.make_event_def('reg1_bd', 1)
game4.setup_gen()
game4.go()
exit_pt = game4.traj.getEvents('exit_ev_reg1_bd')[0]
game2.set(exit_pt[['vx','vy']], exit_pt[['x','y']], by_vel=True)
game2.run()
game2.graphics_refresh(cla=False)
print("Mismatch of pericenter prediction with reduction = %.3f" % abs(min(w2.dist_to_1_vs_peri)))
