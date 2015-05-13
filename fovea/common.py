"""
Common utilities for fovea diagnostics.

"""
from __future__ import division, absolute_import

from numpy import tan
import os, sys

import PyDSTool as dst
from PyDSTool.Toolbox.phaseplane import bisection


class gen_versioner(object):
    _targetlangs = {'dopri': 'c', 'radau': 'c',
                    'vode': 'python', 'euler': 'python'}

    def __init__(self, cwd, model_name, name_base, gen_type, gen_version=0):
        """
        Internal utility to manage versions of Generator objects within single
        session.

        cwd = current working directory (string)

        Option to set known gen version # to reuse:
        Version 0 means not yet created.
        Works across saved and restarted sessions
        """
        self.cwd = cwd
        self.model_name = model_name
        self.name_base = name_base
        self.gen_version = gen_version
        self.gen_type = gen_type
        # keyed by version
        self.used_dsargs = {}
        self.logfile = os.path.join(self.cwd, 'models', 'gen_dsargs_log.sav')
        if os.path.exists(self.logfile):
            # reload previous list of versions
            self.used_dsargs = dst.loadObjects(self.logfile)[0]
        self.classes = {'vode': dst.Generator.Vode_ODEsystem,
                        'dopri': dst.Generator.Dopri_ODEsystem,
                        'radau': dst.Generator.Radau_ODEsystem,
                        'euler': dst.Generator.Euler_ODEsystem}
        self.targetlang = self._targetlangs[gen_type]

    def make(self, dsargs):
        # use cache if available
        for gen_ver, prev_dsargs in self.used_dsargs.items():
            if dst.filteredDict(dsargs, 'name', True) == dst.filteredDict(prev_dsargs,
                                                                  'name', True):
                # compare everything but the name, but check all up to final '_ver<X>'
                parts1 = dsargs.name.split('_')
                parts2 = prev_dsargs.name.split('_') # will have one more part
                if parts1 == parts2[:-2]:
                    print("Reusing identical build")
                    return dst.loadObjects(os.path.join(self.cwd, 'models',
                                                    prev_dsargs.name+'.sav'))[0]
        # no matches
        return self.build(dsargs)

    def build(self, dsargs, is_stiff=False):
        # re-compute in case gen type has been changed
        self.targetlang = self._targetlangs[self.gen_type]
        if is_stiff and self.targetlang == 'python' and self.gen_type == 'vode':
            dsargs.algparams['stiff'] = True

        name = dsargs.name
        if self.gen_version == 0:
            self.gen_version = 1

        # assume it's sufficient to check if .sav file there rather than .so
        found_new = False
        while not found_new:
            filename = os.path.join(self.cwd, 'models', name + '_' + \
                                    self.gen_type + \
                                    '_ver%i'%self.gen_version+'.sav')
            if not os.path.exists(filename):
                found_new = True
            else:
                print(filename + ' already exists')
                self.gen_version += 1
        dsargs.name = name+'_'+self.gen_type+'_ver%i'%self.gen_version
        gen = self.classes[self.gen_type](dsargs)
        model = dst.embed(gen, name=self.model_name, dsi_name='gen')
        self.used_dsargs[self.gen_version] = dsargs.copy()
        self.save_gen(model, name)
        return model

    def load_gen(self, name):
        if self.gen_version == 0:
            raise ValueError("No current version known: set gen_version")
        return dst.loadObjects(os.path.join(self.cwd, 'models', name + '_' + \
                                            self.gen_type + \
                                            '_ver%i'%self.gen_version+'.sav'))[0]

    def save_gen(self, model, name):
        dst.saveObjects(model, os.path.join(self.cwd, 'models', name + '_' + \
                                            self.gen_type + \
                                            '_ver%i'%self.gen_version+'.sav'))
        dst.saveObjects(self.used_dsargs, self.logfile, force=True)


def fdiff(pts, ix):
    """Simple forward finite difference utility:
    Calculates approx to derivative of all coordinates in a pointset at index ix,
    using the time values in the pointset (returning a Point)
    """
    a = pts[ix]
    b = pts[ix+1]
    ts = pts.indepvararray
    ta = ts[ix]
    tb = ts[ix+1]
    return (b-a)/(tb-ta)


class queue(object):
    """
    Create a queue of items that has a limited size.
    As more items than the limit are added, oldest items are
    removed from the end of the queue.
    """

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
    """
    Class to create a queue that has a lookup table.
    """

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
        """Appended item must be a mapping object (dict-like)
        """
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

def checked_scale(sc):
    """Internal utility to verify linear scale syntax:
    None or [None, (ylo,yhi)] or [(xlo,xhi), None]
    of [(xlo,xhi), (ylo,yhi)]
    """
    if sc is not None:
        if len(sc) != 2:
            raise ValueError("Invalid argument for axis scales")
        if sc[0] is not None:
            if len(sc[0]) != 2:
                raise ValueError("X-axis scale must have two components")
        if sc[1] is not None:
            if len(sc[1]) != 2:
                raise ValueError("Y-axis scale must have two components")
    return sc


def castNullArray(null):
    """
    Convert nullcline array data to a list of pairs of arrays usable by plotter.
    """
    return [null[:,0], null[:,1]]

def castNull(null):
    """
    Convert nullcline object to a set of points usable by plotter.
    """
    return [null.array[:,0], null.array[:,1]]
