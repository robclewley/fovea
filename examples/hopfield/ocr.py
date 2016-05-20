from hopfield_network import *

def to_pattern(letter):
    from numpy import array
    return array([+1 if c=='X' else -1 for c in letter.replace('\n','')])

def display(pattern):
    from pylab import imshow, cm, show
    imshow(pattern.reshape((5,5)),cmap=cm.binary, interpolation='nearest')
    show()

A = """
.XXX.
X...X
XXXXX
X...X
X...X
"""
 
Z = """
XXXXX
...X.
..X..
.X...
XXXXX
"""

messy_A = """
XXX..
X...X
.XXXX
X...X
.X..X
"""

distorted_z = """
.....
XXXXX
...X.
..X..
.X...
"""
myNet = HopfieldNetwork(25)
display(to_pattern(A))
myNet.train([to_pattern(A), to_pattern(Z)])
display(myNet.learn([to_pattern(messy_A)]))
