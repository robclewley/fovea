"""
Some basic interface tests
"""

from PyDSTool import *
from bombardier import *


density1 = [0.5, 0.5, 0.5, 0.5, 0.0]
radii1 = [0.04, 0.04, 0.03, 0.03, 0.08]
pos1 = [(-0.2, 0.45), (0.2, 0.7), (0.3, 0.4), (-0.14,0.55),
        (0,0)]

game1 = GUIrocket(radii1, density1, pos1, "Demo-game1",
                  axisbgcol='white')
game1.set( (-60, 0.8) )


def circular_orbit_test():
    # around body 3
    # set max t lower for periodic orbit
    game1.gen.set(tdata=[0,2])
    game1.gen.set(algparams={'max_pts': 30000})
    game1.gen.set(algparams={'max_step': 0.001})
    game1.set_planet_data(1, {'d': 0})
    game1.set_planet_data(2, {'d': 0})
    game1.set_planet_data(4, {'d': 0})
    game1.set( (0, sqrt(game1.masses[2]*G/(0.3-0.25))),
               {'x': 0.25, 'y':0.4})
    game1.go()

def line_test():
    test_line = line_GUI(game1, -0.15, 0.8, 0.1, 0.8)
    game1.selected_object.make_event_def('test', 1)
    game1.setup_gen()
    game1.go()


#line_test()
#circular_orbit_test()