"""
Tests for prototype code to track variable values with callbacks
"""
from __future__ import division
import PyDSTool as dst
from PyDSTool.Toolbox import phaseplane as pp
import numpy as np
from matplotlib import pyplot as plt

import fovea
import fovea.graphics as gx
from fovea.graphics import tracker

class LogWrappedFunction(object):
    def __init__(self, function):
        import inspect
        self.function = function
        self.args, self.varargs, self.keywords, self.defaults = inspect.getargspec(function)
        # ArgSpec(args=['gen', 'subdomain', 'n', 'maxsearch', 'eps', 't', 'jac'], varargs=None, keywords=None, defaults=(None, 5, 1000, 1e-08, 0, None))


    def logAndCall(self, *arguments, **namedArguments):
        print("Calling %s with arguments %s and named arguments %s" %\
                      (self.function.func_name, arguments, namedArguments))
        self.function.__call__(*arguments, **namedArguments)

def logwrap(function):
    return LogWrappedFunction(function).logAndCall

@logwrap
def doSomething(spam, eggs, foo, bar):
    print("Doing something totally awesome with %s and %s." % (spam, eggs))


doSomething("beans","rice", foo="wiggity", bar="wack")

# ============================

def track_attribute(calc_con, attr_name):
    """
    Create decorator to track named attribute used in a calculation
    """
    def decorator(fn):
        #obj =
        #tracker.track_list = [getattr(obj, attr_name)]
        calc_con.workspace
        return fn
    return decorator


def test_func_noattr(x, eps=1e-8):
    """
    mock function that would use a tolerance eps and return a numerical
    object that doesn't contain reference to that tolerance
    """
    return x

# this doesn't let us get at the defaults unless we re-specify a default value
# in the wrapper (see logger above)

def wrap_test_func_noattr(x, eps=1e-8):
    x = test_func_noattr(x, eps)
    res = dst.args(x=x, eps=eps)
    return res

#@track_attribute(cc, 'eps')
def test_func_attr(x, eps=1e-8):
    """
    mock function that would use a tolerance eps and return a numerical
    object that does contain reference to that tolerance
    """
    res = dst.args(val=x, eps=eps)
    return res

x1 = wrap_test_func_noattr(1.0, 1e-8)
x2 = test_func_attr(3.1, 1e-5)

cc = fovea.calc_context(dst.args(tracked_objects=[],
                                 name='saddles'), 'saddles')  # None would be sim object

wksp = cc.workspace
wksp.x1 = x1
wksp.x2 = x2

tracker(cc, 1, text_metadata='eps')

tracker.show()