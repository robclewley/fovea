"""
Common utilities for fovea diagnostics.

"""

from numpy import tan
import os, sys

from PyDSTool.Toolbox.phaseplane import bisection


def progressBar(i, total, width=50):
    percent = float(i)/total
    dots = int(percent*width)
    #os.system('cls' if os.name=='nt' else 'clear')
    if percent > 0:
        sys.stdout.write('\b'*total)
    progress = str('[').ljust(dots+1, '-')
    #print progress.ljust(width, ' ')+str('] %.2f%%' % (percent*100.)),
    sys.stdout.write(progress.ljust(width, ' ')+str('] %.2f%%' % (percent*100.)))
    sys.stdout.flush()


def fdiff(pts, ix):
    """Simple forward finite difference utility:
    Calculates approx to derivative of all coordinates in a pointset at index ix,
    using the time values in the pointset (returning a Point)"""
    a = pts[ix]
    b = pts[ix+1]
    ts = pts.indepvararray
    ta = ts[ix]
    tb = ts[ix+1]
    return (b-a)/(tb-ta)


def make_mline(pt, angle):
    """Utility factory function to make m a function of V (constrained along a line
    starting at pt at the given angle): Returns a callable m(V)"""
    def m(V):
        return (V-pt['V'])*tan(angle) + pt['Na.m']
    return m



class queue(object):
    """Create a list of items with a limited size.  As more items are added, older
    items are removed from the end of the stack."""


    def __init__(self, size, data=None):

        self.size = size

        if data is None:
            self.data = []
        elif type(data) is not list:
            raise TypeError('queue object must be a list')
        else:
            assert len(data) <= size, 'Initial data must be smaller than queue size'
            self.data = data

    def __call__(self, i=None):
        if i is None:
            ret = self.data
        else:
            ret = self.data[i]
        return ret

    def append(self, item):

        if len(self.data) < self.size:
            print len(self.data)
            self.data.append(item)
            return []
        else:
            self.data.append(item)
            extra = len(self.data) - self.size
            popped = []
            for i in range(extra):
                popped.append(self.data.pop(i))
            return popped


class hashedQ(object):
    """ Class to create a queue that has a lookup table."""

    def __init__(self, size):
        self.size = size
        self.keys = []
        self.data = []

    def __call__(self, key=None):
        if key is None:
            ret = {}
            for key in self.keys:
                ix = self.keys.index(key)
                ret.update({key: self.data[ix]})
        else:
            ix = self.keys.index(key)
            ret = self.data[ix]
        return ret

    def append(self, item):

        assert type(item) is dict, "Appended item must be of type dict"

        for key in item:
            self.keys.append(key)
            self.data.append(item[key])

        if len(self.keys) > self.size:
            extra = len(self.keys) - self.size
            for i in range(extra):
                self.keys.pop(i)
                self.data.pop(i)
            return False

        return True


def castNullArray(null):
    """
    Convert nullcline array data to a list of pairs of arrays useable by plotter.
    """
    return [null[:,0], null[:,1]]

def castNull(null):
    """
    Convert nullcline object to a set of points useable by plotter.
    """
    return [null.array[:,0], null.array[:,1]]
