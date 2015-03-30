"""
Some solutions to hit target opposite

"""
from PyDSTool import *
from bombardier import *


# ================================

sol1 = (20.22, 0.646) # loop each side
sol1b = (17.58, 0.808) # loop each side
sol2 = (10.37, 0.8) # no loops
sol3 = (17.58, 0.773) # right-sided loop
sol4 = (1.934, 0.528) # left-sided loop
sol5 = (4.57, 0.38) # left-sided loop

end_sol1 = (4.3835, 1.884)
end_sol1_ic = {'x': -0.03564, 'y': 0.55952}

game1.set(sol1)
game1.go()

ptdata = findpt(game1, 0.4, .25, .6)
pt = ptdata[2]
game1.set( (pt['bearing'], pt['speed']),
           pt['x','y'])
#game1.set(end_sol1, end_sol1_ic)
#game1.setAng(4.868)
#game1.setVel(0.5458)

# balance line through first quartile of sol1
bl = line_GUI(game1, -0.2,0.5,0.4,0.4)
# need to ensure (x,y) on bl
print project(game1, 0.1, 0.45, 3, bl)
