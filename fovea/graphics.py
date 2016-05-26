"""
Graphical user interface and plotting tools for dynamical systems

Rob Clewley, 2015
based on work by:
Bryce Chung and Rob Clewley, 2012


==== Usage notes:
Plotting styles can be given either as a string or as a dictionary of
  style kwargs suitable for plot command

"""
from __future__ import division, absolute_import

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RectangleSelector
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from copy import copy
from math import *
from collections import OrderedDict
import hashlib, time
import euclid as euc

import warnings

from PyDSTool import * #Need Events from here.

from PyDSTool import args, numeric_to_traj, Point, Points
import PyDSTool.Toolbox.phaseplane as pp
# for potentially generalizable functions and classes to use
import PyDSTool as dst

# local imports
from fovea.common import *
import fovea.domain2D as dom
from fovea.diagnostics import get_unique_name

# ----------------------------------------------

#Test python version
try:
    input = raw_input
except NameError:
    # python 3
    pass


def force_line_to_extent(a, b, p_domain, coordnames):
    """
    The infinite line through points *a* and *b*
    cuts the rectangular domain at the two points returned.

    a, b are Point2D objects for the line endpoints.
    """
    x, y = coordnames
    # line is a euclid.Line2 (infinite in both directions)
    a_pt = euc.Point2(a[x],a[y])
    b_pt = euc.Point2(b[x],b[y])
    #plt.plot(np.array([a_pt, b_pt]).T[0],
    #         np.array([a_pt, b_pt]).T[1], 'r')
    line = euc.Line2(a_pt, b_pt)
    # test each of four sides as finite line segments
    combos = ( (0,0,0,1), (0,0,1,0), (1,1,0,1), (1,1,1,0) )
    p1 = p2 = None
    for xi1, yi1, xi2, yi2 in combos:
        bd1 = euc.Point2(p_domain[x][xi1], p_domain[y][yi1])
        bd2 = euc.Point2(p_domain[x][xi2], p_domain[y][yi2])
        #plt.plot(np.array([bd1, bd2]).T[0],
        #         np.array([bd1, bd2]).T[1], 'k')
        border_seg = euc.LineSegment2(bd1, bd2)
        pt = line.intersect(border_seg)
        if pt is not None:
            # intersection exists
            # truncate or extend whichever point of line is closer to
            # border segment
            dist1 = border_seg.distance(line.p1)
            dist2 = border_seg.distance(line.p2)
            if dist1 > dist2:
                if p2 is None:
                    p2 = pt
                else:
                    p1 = pt
            else:
                if p1 is None:
                    p1 = pt
                else:
                    p2 = pt
    if p1 is not None and p2 is not None:
        return pp.Point2D(p1), pp.Point2D(p2)
    else:
        raise ValueError("No intersection")


class Plotter(object):

    colors = ['b', 'g', 'r', 'c', 'm', 'k', 'y']

    def __init__(self, dm=None):
        """
        Optionally pass in a diagnostic manager object as `dm`.
        """
        self.clean()
        self.dm = dm
        self.save_status = False
        self.wait_status = False

    def clean(self):
        """
        Delete all figure data (doesn't clear figure displays)
        """
        self.figs = {}
        self._max_fig_num = 0
        self.active_layer = None
        self.active_layer_structs = None
        # lookup dictionary mapping axes handles to layer name lists
        # and sub-plot index strings
        self.subplot_lookup = {}
        self.currFig = None
        # record whether this class ever called show()
        self.shown = False

    def auto_scale_domain(self, xcushion=0, ycushion=0, subplot=None, figure=None):
        """
        Set domain limits to that of the data in all layers
        with the greatest extent.

        cushion is the fraction of the total x/y domain to extend the domain (leaving some white space).
        """
        # ISSUE: Setting domains at the figure level (as opposed to subplot level) doesn't work.
        # Changes will likely need to be made in other functions, not here (such as build_layers or build_plotter)
        x_extent = [0, 0]
        y_extent = [0, 0]

        def auto_minmax(x, f):
            """
            Scalar inputs are automatically the min/max.
            """
            try:
                out = f(x)
            except TypeError:
                out = x
            return out

        if not figure: 
            figure = self.currFig

        found_fig = False
        for figName, fig in self.figs.items():
            if figure == figName:
                found_fig = True
            # Extract only those layers in the chosen subplot.
            if subplot:
                subplot_struct = fig.arrange[subplot]
                subplot_layers = subplot_struct['layers']
                layer_info = { layer_name : fig.layers[layer_name] for layer_name in subplot_layers }
            else:
                layer_info = fig.layers
            for layerName, layer in layer_info.items():
                if layer.kind != 'text':
                   for dName, d in list(layer['data'].items()):
                       data_points = d['data']
                       x_extent[0] = min(auto_minmax(data_points[0], min), x_extent[0])
                       x_extent[1] = max(auto_minmax(data_points[0], max), x_extent[1])
                       y_extent[0] = min(auto_minmax(data_points[1], min), y_extent[0])
                       y_extent[1] = max(auto_minmax(data_points[1], max), y_extent[1])

        if not found_fig:
            raise ValueError("No such figure")

        if subplot:
            x_length = x_extent[1] - x_extent[0]
            y_length = y_extent[1] - y_extent[0]
            subplot_struct['scale'] = [tuple([x_extent[0] - xcushion*x_length, x_extent[1] + xcushion*x_length]),
                                             tuple([y_extent[0] - ycushion*y_length, y_extent[1] + ycushion*y_length])]
        else:
            fig.domain = (x_extent, y_extent)
            plt.figure(fig.fignum)
            plt.xlim(x_extent)
            plt.ylim(y_extent)

    def set_active_layer(self, layer, figure=None):
        """
        Sets the active_layer attribute to be a pair
        of the given figure name (optional, defaults to Master)
        and the named layer struct.
        """
        fig_struct, figure = self._resolve_fig(figure)
        try:
            layer_struct = fig_struct.layers[layer]
        except KeyError:
            raise KeyError("Unknown layer: %s" % layer)
        else:
            self.active_layer = (figure, layer)
            self.active_layer_structs = (fig_struct, layer_struct)

    def show_legends(self, figure=None, subplot=None):
        """
        Show all figure legends of visible data layers to stdout.
        Option to filter results by sub-plot string code or name, e.g. '21'
        """
        fig_struct, figure = self._resolve_fig(figure)
        arrPlots = fig_struct.arrange

        if not arrPlots:
            print("  This figure has no associated layer data.")
            return
        for ixstr, spec in arrPlots.items():
            if subplot is not None and subplot not in \
                              [ixstr, spec['name']]:
                continue
            layer_info = spec['layers']
            for layer in layer_info:
                lay = fig_struct.layers[layer] #Doesn't interpret layer = * properly
                if not lay.display or lay.kind != 'data':
                    continue
                print("\nLAYER: %s" % str(layer))
                print("  sub-plot: %s" % ixstr)
                print("  style %s" % lay.style)
                if lay.axes_vars:
                    print("  axes: %s - %s" % (lay.axes_vars[0], lay.axes_vars[1]))
                print("  data:")
                for dname, dstruct in lay.data.items():
                    if not dstruct['display']:
                        continue
                    print("     name: %s, style: %s" % (dname, dstruct['style']))

    ## Figure Management Tools ##

    def add_fig(self, label, title="", xlabel="", ylabel="", tdom=None,
               domain=None, display=True):
        """
        User can add figures to plotter for data shown in different windows

        label        String to specify name of figure
        domain       Pair of axis domain extent value pairs
        title        String to specify title to display on figure
        xlabel       String to specify label on x-axis
        ylabel       String to specify label on y-axis
        tdom         time domain (if applicable, and optional)
        display      Setting to determine whether or not figure is plotted when
                       plotter is built (defaults to True)
        """
        # tdom is possibly redundant (or too specialized?)
        if domain is None:
            raise ValueError("Must supply domain for figure")

        # Check to see if label already exists
        if label in self.figs:
            raise KeyError("Figure label already exists!")

        self.currFig = label

        figAttr = args()
        figAttr.title = title
        figAttr.xlabel = xlabel
        figAttr.ylabel = ylabel
        figAttr.display = display
        self._max_fig_num += 1
        figAttr.fignum = self._max_fig_num
        figAttr.layers = {}
        figAttr.shape = [1,1]
        figAttr.arrange = []
        figAttr.window = None
        figAttr.autoscaling = True
        figAttr.tdom = tdom
        figAttr.domain = domain

        self.figs.update({label: figAttr})


    def copy_fig(self, newFig, oldFig):
        """
        Duplicates all figure, layer, and data information without arrangement
        information.
        """
        # Check to see if label already exists #
        if newFig in self.figs:
            raise KeyError("New figure label already exists!")

        self.currFig = newFig

        old = self.figs[oldFig]

        figAttr = args()
        figAttr.title = old.title
        figAttr.xlabel = old.xlabel
        figAttr.ylabel = old.ylabel
        figAttr.display = old.display
        self._max_fig_num += 1
        figAttr.fignum = self._max_fig_num
        figAttr.shape = [1,1]
        figAttr.arrange = []
        figAttr.window = None
        figAttr.autoscaling = True

        figAttr.layers = {}

        self.figs.update({newFig: figAttr})

        for layer in old.layers:
            oldLay = self.figs[oldFig].layers[layer]
            self.add_layer(layer, display=oldLay.display, zindex=oldLay.zindex, style=oldLay.style)


    def set_fig(self, label=None, **kwargs):
        """
        Sets current figure given by parameter 'label', if given.

        Also alows user to set properties of the named (or current) figure.
        See add_fig() for properties to set
        """
        # Check to see that figure exists #
        if label == None and self.currFig != None:
            label = self.currFig
        elif label == None and self.currFig == None:
            raise ValueError("Must set current figure")
        elif label != None:
            if label not in self.figs:
                raise KeyError("Figure does not exist with specified label!")

        self.currFig = label

        for kw in kwargs:
            # Check to see that keyword is valid #
            if kw not in self.figs[label]:
                raise KeyError("Unable to set figure property: Invalid keyword argument")

            self.figs[label][kw] = kwargs[kw]


    def clear_fig_data(self, figure_name):
        """
        Clear all figure data for named figure
        Ignore if figure name is found.
        """
        if figure_name in self.figs:
            for layer in self.figs[figure_name].layers:
                self.figs[figure_name].layers[layer].data = {}

            if figure_name == self.active_layer[0]:
                self.active_layer = None
                self.active_layer_structs = None

            if self.currFig == figure_name:
                self.currFig = None


    ## Layer Management Tools ##

    def arrange_fig(self, shape, arrPlots, figure=None):
        """
        Separate layers in a figure into different subplots.

        shape      [rows,cols] where rows is number of rows starting from 1
                   and cols is number of cols starting from 1
        arrPlots   Dict of dicts of subplots indexed by row-col.  Each
                   subplot dict has the following keys:

                   name        Name of the subplot to appear in figure
                   layers      List of layers to be displayed in subplot
                   axes_vars   Names of axes, [xlabel, ylabel]
                   scale       Scale of axes [(xlo, xhi), (ylo, yhi)]
                                (will NOT overwrite any already declared for a layer)

        figure     Figure name to apply changes to. If left blank, current
                   figure is used

        """
        fig_struct, figure = self._resolve_fig(figure)

        if len(shape) != 2:
            raise ValueError("shape must be (rows,cols)")

        # make copy of arrPlots in case change singleton layer names to lists
        arrPlots = copy(arrPlots)

        #Ensure subplot positions are consistent with figure shape.
        for ixstr, spec in arrPlots.items():
            if int(ixstr[0])*int(ixstr[1]) > shape[0]*shape[1]:
                raise ValueError("Position does not exist in subplot arrangement.")

            layer_info = spec['layers']
            if len(spec['axes_vars']) > 3:
                raise ValueError("Cannot have more than three axis titles.")
            else:
                axes_vars = spec['axes_vars']

        fig_struct.shape = shape
        fig_struct.arrange = arrPlots


    def add_layer(self, layer_name, figure=None, set_to_active=True, subplot=None,
                 **kwargs):
        """
        User method to add data sets to a layer in a figure.
        By default, the latest added layer becomes the 'active' layer,
        unless set_to_active=False.

        figure  name
        layer   name
        display toggle Boolean
        kind    a user-defined kind, e.g. 'data', 'vline', 'hline', 'epoch',
                     'text', etc.
        scale   a pair of axes scale pairs or None
        zindex  (currently unused)
        """
        # ISSUE: Not sure that figure or whole layer display attribute values
        # are respected
        fig_struct, figure = self._resolve_fig(figure)

        # Check to see layer does not already exist
        if layer_name in fig_struct.layers:
            raise KeyError("Layer name already exists in figure!")

        # default plot style generated from list of colors
        color = self.colors[len(fig_struct.layers)%len(self.colors)]
        line = '-'
        style = color+line

        layAttrs = args()
        layAttrs.data = {}
        layAttrs.display = True
        layAttrs.force = False
        # ISSUE: zindex ordering not actually used
        # default zindex based on order of layer declaration
        layAttrs.zindex = len(fig_struct.layers)+1
        layAttrs.style = style
        #layAttrs.axes = None
        layAttrs.scale = None
        layAttrs.kind = 'data'
        layAttrs.dynamic = False
        layAttrs.trajs = {}
        layAttrs.axes_vars = []
        layAttrs.handles = OrderedDict({})
        #layAttrs.linewidth = None

        for kw in kwargs:
            # Check to see that parameter exists in layers
            # ISSUE: Possibly change to account for different properties of
            # specific artist objects?
            if kw not in layAttrs:
                raise KeyError("Parameter is not a property of the layer.")
            layAttrs[kw] = kwargs[kw]

        fig_struct.layers[layer_name] = layAttrs
        if set_to_active:
            self.set_active_layer(layer_name, figure=figure)

        #Adds the layer to the arrange, translating an axes 'subplot' into a subplot string.
        if subplot is None:
            pass
        elif isinstance(subplot, str):
            fig_struct['arrange'][subplot]['layers'] += [layer_name]
        elif isinstance(subplot, mpl.axes.Axes):
            for key in fig_struct['arrange'].keys():
                if subplot is fig_struct['arrange'][key]['axes_obj']:
                    fig_struct['arrange'][key]['layers'] += [layer_name]
                    break
        else:
            raise TypeError("subplot must be either string or axes object.")


    def set_layer(self, label, figure=None, **kwargs):
        """
        Arrange data sets in a figure's layer
        """
        # figure will be the same before and after unless figure was
        # None, in which case defaults to name of master figure
        fig_struct, figure = self._resolve_fig(figure)

        # Check to see that layer exists
        if label not in fig_struct.layers:
            raise KeyError("Layer does not exist in figure!")

        for kw in kwargs:
            # Check to see that parameter exists in layers
            # ISSUE: Possibly change to account for different properties of
            # specific artist objects?
            if kw not in fig_struct.layers[label]:
                raise KeyError("Parameter is not a property of the layer.")

            fig_struct.layers[label][kw] = kwargs[kw]

    def add_patch(self, data, patch, figure=None, layer=None, subplot=None, name=None,
                display=True, force=False, log=None, **kwargs):
        """
        patch is a matplotlib patch class. Accepts kwargs for patch objects.
        """
        try:
            size = np.shape(data)
        except:
            raise TypeError("Data must be castable to a numpy array")

        # Check to see that there is an x- and y- dataset
        try:
            if not size[0] == 2:
                raise ValueError("Data must contain 2 or 3 seqs of data points")
        except IndexError:
            pass

        fig_struct, figure = self._resolve_fig(figure)
        if layer is None:
            layer = self.active_layer[1]
            layer_struct = self.active_layer_structs[1]
        else:
            layer_struct = self._resolve_layer(figure, layer)

        if not layer_struct.kind == 'patch':
            raise ValueError("Incompatible layer type (should be `patch`)")

        # d is a dictionary mapping 'name' to a dictionary of fields for the
        # numerical plot data (key 'data'), style (key 'style'), and display
        # boolean (key 'display').
        #
        # The numerical data for d is given by the method argument also called data
        d = layer_struct.data

        # Check to see if data name already exists
        if name in d and not force:
            raise KeyError("Data name already exists in layer: %s" %name)

        # Create name if none is provided
        if name is None:
            name = get_unique_name(figure+'_'+layer)

        if log:
            log.msg("Added plot data", figure=figure, layer=layer, name=name)

        kwargs.update({'data':data,'patch':patch, 'display':display, 'subplot':subplot})
        d.update({name: kwargs})
        layer_struct.force = force

    def add_obj(self, data, obj, figure=None, layer=None, subplot=None, name=None,
                 display=True, force=False, log=None, linewidth=1, markersize=6, **kwargs):
        # ISSUE: May want to combine this with patch, e.g. to addClass,
        # or keep this as add_context_obj for an internal function (not for end user).
        try:
            size = np.shape(data)
        except:
            raise TypeError("Data must be castable to a numpy array")

        # Check to see that there is an x- and y- dataset
        try:
            if not size[0] == 2:
                raise ValueError("Data must contain 2 or 3 seqs of data points")
        except IndexError:
            pass

        fig_struct, figure = self._resolve_fig(figure)
        if layer is None:
            layer = self.active_layer[1]
            layer_struct = self.active_layer_structs[1]
        else:
            layer_struct = self._resolve_layer(figure, layer)

        if not layer_struct.kind == 'obj':
            raise ValueError("Incompatible layer type (should be `obj`)")

        # d is a dictionary mapping 'name' to a dictionary of fields for the
        # numerical plot data (key 'data'), style (key 'style'), and display
        # boolean (key 'display').
        #
        # The numerical data for d is given by the method argument also called data
        d = layer_struct.data

        # Check to see if data name already exists
        if name in d and not force:
            raise KeyError("Data name already exists in layer: %s" %name)

        # Create name if none is provided
        if name is None:
            name = get_unique_name(figure+'_'+layer)

        if log:
            log.msg("Added plot data", figure=figure, layer=layer, name=name)

        #Translate axes into subplot string.
        if isinstance(subplot, mpl.axes.Axes):
            for key in fig_struct['arrange'].keys():
                #Direct comparison of objects was unsuccessful. Perhaps caused by rebuilding plotter. Comparing titles works.
                if subplot.get_title() is fig_struct['arrange'][key]['axes_obj'].get_title():
                    subplot = key
                    break

        kwargs.update({'data':data,'obj':obj, 'display':display, 'subplot':subplot,
                       'linewidth':linewidth, 'markersize':markersize, 'selected':False})
        d.update({name: kwargs})
        layer_struct.force = force

    def add_data(self, data, figure=None, layer=None, subplot=None, style=None, linewidth = 1,
                markersize= 6, zorder= 1, name=None, display=True, force=False, traj = None, log=None):
        """
        User tool to add data to a named layer (defaults to current active layer).
        *data* consists of a pair of sequences of x, y data values, in the same
        format as would be passed to matplotlib's plot.

        Use *force* option only if known that existing data must be overwritten.
        Add a diagnostic manager's log attribute to the optional *log* argument
        to have the figure, layer, and data name recorded in the log.

        *display* option (default True) controls whether the data will be
        visible by default.
        """
        # Check to see that data is a list or array
        try:
            size = np.shape(data)
        except:
            raise TypeError("Data must be castable to a numpy array")

        # Check to see that there is an x- and y- (or z-) dataset
        try:
            if size[0] not in (2,3):
                raise ValueError("Data must contain 2 or 3 seqs of data points")
        except IndexError:
            pass

        fig_struct, figure = self._resolve_fig(figure)
        if layer is None:
            layer = self.active_layer[1]
            layer_struct = self.active_layer_structs[1]
        else:
            layer_struct = self._resolve_layer(figure, layer)

        # inherit default style from layer if not given
        if style is None:
            style = layer_struct.style

        # d is a dictionary mapping 'name' to a dictionary of fields for the
        # numerical plot data (key 'data'), style (key 'style'), and display
        # boolean (key 'display').
        #
        # The numerical data for d is given by the method argument also called data
        d = layer_struct.data

        # Check to see if data name already exists
        if name in d and not force:
            raise KeyError("Data name already exists in layer: %s" %name)

        # Create name if none is provided
        if name is None:
            name = get_unique_name(figure+'_'+layer)

        if log:
            log.msg("Added plot data", figure=figure, layer=layer, name=name)
        d.update({name: {'data': data, 'style': style, 'linewidth':linewidth, 'markersize':markersize,
                         'zorder':zorder,'display': display, 'subplot': subplot, 'selected':False}})
        layer_struct.force = force

        # ISSUE: _update_traj only meaningful for time-param'd trajectories
        # Maybe a different, more general purpose solution is needed

        self._update_traj(figure, layer, traj = traj)

    def set_point(self, name, pt, layer, figure=None):
        """
        Set point coordinates in given layer, specified by pt.

        Figure defaults to currFig if not specified.
        """
        fig_struct, figure = self._resolve_fig(figure)
        lay = fig_struct.layers[layer]
        pt_struct = lay.data[name]
        try:
            lay.handles[name].set_data(pt)
        except KeyError:
            # object not actually plotted yet
            pass

    def set_data_2(self, label, layer, figure=None, **kwargs):
        """
        Set properties for a specific dataset in given layer, specified
        by the keys of the keyword arguments

        ISSUE: Original set_data doesn't seem to provide any additional functionality beyond what set_layer already does.
        This function is meant to behave as set_layer but at the 'data' level in Fovea's object hierarchy.
        """
        # figure will be the same before and after unless figure was
        # None, in which case defaults to name of master figure
        fig_struct, figure = self._resolve_fig(figure)

        # Check to see that layer exists
        if layer not in fig_struct.layers:
            raise KeyError("Layer does not exist in figure!")

        # Check to see that data exists
        if label not in fig_struct.layers[layer].data:
            raise KeyError("Data does not exist in layer!")

        for kw in kwargs:
            # Check to see that parameter exists in layers
            # ISSUE: Possibly change to account for different properties of
            # specific artist objects?
            if kw not in fig_struct.layers[layer].data[label]:
                raise KeyError("Parameter is not a property of the data.")

            fig_struct.layers[layer].data[label][kw] = kwargs[kw]


    def set_data(self, layer, figure=None, **kwargs):
        """
        Set data properties in given layer, specified by
        the keys of the keyword arguments.

        Figure defaults to currFig if not specified.
        """
        fig_struct, figure = self._resolve_fig(figure)

        # Check to see that layer exists
        try:
            layer_struct = fig_struct.layers[layer]
        except KeyError:
            raise KeyError("Layer does not exist in figure!")

        for kw in kwargs:
            # ISSUE: Possibly change to account for different properties
            # of specific artist objects?
            try:
                # for dict updates, e.g. data={<line_name>: {'data': [<array_data>]}}
                layer_struct[kw].update(kwargs[kw]) # = kwargs[kw]
            except AttributeError:
                # for non-dict updates, e.g. display=True
                layer_struct[kw] = kwargs[kw]
            except KeyError:
                raise KeyError("Parameter '%s' is not a property of the data." %kw)
            else:
                # dict update: check whether data was changed and propogate to plot
                # object
                if 'data' in kwargs:
                    for objname, objdata in kwargs['data'].items():
                        try:
                            obj = layer_struct.handles[objname]
                        except KeyError:
                            # object not actually plotted yet
                            pass
                        #else:
                            #print("OBJ: ", obj)
                            #obj.set_data(objdata['data'])

        self._update_traj(figure, layer)



    def _update_traj(self, figure_name, layer_name, traj= None):
        """
        Create an interpolated trajectory from the data in the layer
        This may no longer be necessary (it's not general purpose for fovea)
        """
        fig_struct = self.figs[figure_name]
        layer_struct = fig_struct.layers[layer_name]

        #pointset_to_traj requires an independent variable.
        try:
            if traj.indepvartype is None:
                traj = None
        except AttributeError:
            pass

        ##ISSUE: Should this loop through all the dstructs? Seems like it just needs the one associated
        ## with the call to add_data.
        for name, dstruct in layer_struct.data.items():
            if traj is not None:
                layer_struct.trajs[name] = pointset_to_traj(traj)
                layer_struct.trajs[name].name = layer_name+'.'+name
            # catch for when the data is individual points
            else:
                try:
                    layer_struct.trajs[name] = numeric_to_traj([dstruct['data'][0], dstruct['data'][1]],
                                            layer_name+'.'+name,
                                            coordnames=['x', 'y'],
                                            indepvar=dstruct['data'][0],
                                            discrete=False)
                except ValueError:
                    pass


    def toggle_display(self, names=None, layer=None, figure=None, log=None):
        """
        Toggle the display attribute of an object or list of objects, a whole
        layer, or a whole figure, depending on which optional arguments are
        left as None.

        If figure is left None (default), then the current figure is used.

        Object *names* can be a singleton string or a list of strings in the
        same layer. If these names are provided, layer can be left None to
        select default layer.
        """
        fig_struct, figure = self._resolve_fig(figure)

        if names is None:
            if layer is None:
                # toggle whole figure display
                fig_struct.display = not fig_struct.display
            else:
                layer_struct = self._resolve_layer(figure, layer)
                # toggle whole layer display
                layer_struct.display = not layer_struct.display
        else:
            layer_struct = self._resolve_layer(figure, layer)
            if isinstance(names, str):
                names = [name]
            for name in names:
                disp = lay.data[name]['display']
                lay.data[name]['display'] = not disp

    def set_display(self, names, display, layer, figure=None):
        """
        *names* can be a singleton string or a list of strings in the same
        layer.
        """
        fig_struct, figure = self._resolve_fig(figure)
        try:
            lay = fig_struct.layers[layer]
        except KeyError:
            raise ValueError("Layer %s does not exist in figure %s"%(layer, figure))

        assert display in [True, False]
        if isinstance(names, str):
            names = [name]
        for name in names:
            lay.data[name]['display'] = display

    def append_data(self, data, layer, name, figure=None, log=None):
        """
        Append data to the given named data in the given layer.

        display attribute of existing data will continue to apply.
        ISSUE: Doc string?
        """
        fig_struct, figure = self._resolve_fig(figure)
        try:
            lay = fig_struct.layers[layer]
        except KeyError:
            raise ValueError("Layer %s does not exist in figure %s"%(layer, figure))

        try:
            dataset = lay.data[name]
        except KeyError:
            raise ValueError("Dataset %s does not exist in layer" % name)

        if isinstance(data, pp.Point2D):
            x = data.x
            y = data.y
        else:
            try:
                x = data[0]
                y = data[1]
                assert len(data) == 2
            except:
                raise TypeError("Point must be of type Point2D or iterable")

        dataset['data'][0].append(x)
        dataset['data'][1].append(y)

        if log:
            log.msg("Appended plot data", figure=figure, layer=layer,
                        name=name)

        self._update_traj(layer, figure)


    def add_line_by_points(self, pts, figure=None, layer=None, style=None,
                        name=None, display=True, log=None):
        """
        Add line based on two given Point2D points
        """
        pt1, pt2 = pts
        try:
            self.add_data([[pt1.x, pt2.x], [pt1.y, pt2.y]], figure=figure, layer=layer, style=style, name=name,
                         display=display, log=log)
        except AttributeError:
            raise ValueError("First argument must be [Point2D, Point2D]")


    def add_point(self, pt, figure=None, layer=None, style=None, name=None,
                 display=True, log=None):
        """
        Add single Point2D point. style argument should include specification of
        marker symbol and color, if desired.
        """
        self.add_data([[pt.x], [pt.y]], figure=figure, layer=layer, style=style,
                     name=name, display=display, log=log)


    def add_vline(self, x, figure=None, layer=None, subplot=None, style=None, name='vline',
                 log=None):
        """
        Add vertical line.
        """
        # ISSUE: Same issue as add_hline -- see below
        fig_struct, figure = self._resolve_fig(figure)
        layer_struct = fig_struct.layers[layer]
        sc = layer_struct.scale
        if sc is None:
            ydom = self.figs[figure].domain[1]
        else:
            if sc[1] is None:
                ydom = self.figs[figure].domain[1]
            else:
                ydom = sc[1]
        self.add_data([[x, x], ydom], figure=figure, layer=layer, subplot=subplot,
                         style=style, name=name, log=log)


    def add_hline(self, y, figure=None, layer=None, style=None, name='hline',
                 log=None):
        """
        Add horizontal line.
        """
        # ISSUE: This should be changed to use ax.axhline, which automatically
        # always spans the x axis with the default coords settings.
        fig_struct, figure = self._resolve_fig(figure)
        layer_struct = fig_struct.layers[layer]
        sc = layer_struct.scale
        if sc is None:
            xdom = self.figs[figure].domain[0]
        else:
            if sc[0] is None:
                xdom = self.figs[figure].domain[0]
            else:
                xdom = sc[0]
        self.add_data([xdom, [y, y]], figure=figure, layer=layer,
                     style=style, name=name, log=log)


    def add_text(self, x, y, text, use_axis_coords=False, name=None, layer=None, subplot=None,
                figure=None, display=True, style=None, force=False, log=None):
        """
        Use style to select color (defaults to black).
        """
        # ISSUE: Cannot set font size or weight
        fig_struct, figure = self._resolve_fig(figure)
        if layer is None:
            layer = self.active_layer[1]
            layer_struct = self.active_layer_structs[1]
        else:
            layer_struct = self._resolve_layer(figure, layer)

        if not layer_struct.kind == 'text':
            raise ValueError("Incompatible layer type (should be `text`)")

        # The numerical data for d is given by the method argument also called data
        d = layer_struct.data

        # Check to see if data name already exists
        if name in d and not force:
            raise KeyError("Data name already exists in layer: %s" %name)

        # Create name if none is provided
        if name is None:
            name = get_unique_name(figure+'_'+layer)

        d.update({name: {'data': [x, y], 'text': text, 'display': display,
                         'style': style, 'use_axis_coords': use_axis_coords, 'subplot': subplot}})

        if log:
            log.msg("Added text data", figure=figure, layer=layer, name=name)


    def set_text(self, name, text, layer, pos=None, figure=None):
        """
        Use optional `pos` argument to set new position (coord value pair).
        """
        fig_struct, figure = self._resolve_fig(figure)
        lay = fig_struct.layers[layer]
        text_struct = lay.data[name]
        text_struct['text'] = text
        try:
            lay.handles[name].set_text(text)
        except KeyError:
            # name not actually plotted yet
            pass
        if pos is not None:
            text_struct['data'] = pos


    def _subplots(self, layers, fig_name, rebuild=False):
        fig_struct = self.figs[fig_name]
        fig = plt.figure(fig_struct.fignum)
        if rebuild:
            for axs in fig.get_axes():
                if axs not in [bttn.ax for bttn in self.gui.widgets.values()]:
                    fig.delaxes(axs)

        # Build up each subplot, left to right, top to bottom
        shape = fig_struct.shape
        for i in range(shape[0]):
            for j in range(shape[1]):
                ixstr = str(i+1) + str(j+1)
                try:
                    subplot_struct = fig_struct.arrange[ixstr]
                except (KeyError, TypeError):
                    # type error if arrange is empty list (e.g. for shape=[1,1])
                    continue
                layer_info = subplot_struct['layers']
                if not isinstance(layer_info, list):
                    if layer_info == '*':
                        layer_info = list(fig_struct.layers.keys())
                    else:
                        # singleton string layer name
                        layer_info = [layer_info]

                try:
                    scale = subplot_struct['scale']
                except KeyError:
                    subplot_struct['scale'] = None
                    scale = None

                if rebuild:
                    #Check if projection type has been specified for layer.
                    try:
                        ax = fig.add_subplot(shape[0], shape[1], shape[1]*i + j+1, projection= subplot_struct['projection'])
                    except KeyError:
                        ax = fig.add_subplot(shape[0], shape[1], shape[1]*i + j+1)

                    subplot_struct['axes_obj'] = ax
                    ax.set_title(subplot_struct['name'])
                    axes_vars = subplot_struct['axes_vars']
                    ax.set_xlabel(axes_vars[0])
                    ax.set_ylabel(axes_vars[1])

                    if len(axes_vars) == 3:
                        if subplot_struct['projection'] != '3d':
                            raise ValueError("Cannot have 3 axes variables on a layer where projection is not '3d'")
                        ax.set_zlabel(axes_vars[2])

                else:
                    ax = subplot_struct['axes_obj']

                if 'callbacks' in subplot_struct:
                    if ax not in self.gui.cb_axes:
                        self.gui.cb_axes.append(ax)
                        self.gui.RS_boxes[ax] = RectangleSelector(ax, self.gui.onselect_box, drawtype= 'box')
                        self.gui.RS_lines[ax] = RectangleSelector(ax, self.gui.onselect_line, drawtype='line')
                        self.gui.RS_boxes[ax].set_active(False)
                        self.gui.RS_lines[ax].set_active(False)

                if 'legend' in subplot_struct:
                    handles = subplot_struct['legend']
                    subplot_struct['axes_obj'].legend(handles= handles)

                # refresh this in case layer contents have changed
                self.subplot_lookup[ax] = (fig_name, layer_info, ixstr)

                # ISSUE: these should be built into Plotter's figure domains instead
                if scale is not None:
                    # scale may be [None, None], [None, [ylo, yhi]], etc.
                    try:
                        ax.set_xlim(scale[0])
                    except TypeError:
                        pass
                    try:
                        ax.set_ylim(scale[1])
                    except TypeError:
                        pass

                self.build_layers(layer_info, ax, rebuild=rebuild)


    def show(self, update='current', rebuild=False, force_wait=None, ignore_wait= False):
        """
        Apply all figures' domain limits and update/refresh with any
        latest data added or changed.

        By default, will udpate 'current' (or set to None or 'all') layers
        and figures first, in case of newly added/altered data content.

        Optional rebuild = True argument (default False) to clear whatever
        was selected by `update` and rebuild it from scratch.

        Optional force_wait argument to override current wait_status attribute
        to either pause execution until <RETURN> key pressed or stop pausing.

        Set ignore_wait to True if updates from one call to show() must take effect
        while the wait status of a different call is set to True.
        """
        if update == 'current':
            fig_struct, fig_name = self._resolve_fig(None)
            layers = {fig_name: fig_struct.layers.keys()}
            figures = [fig_name]
        elif update == 'all':
            layers = {}
            figures = []
            for fig_name, fig_struct in self.figs.items():
                figures.append(fig_name)
                layers[fig_name] = list(fig_struct.layers.keys())
        if update is not None:
            for fig_name in figures:
                self._subplots(layers[fig_name], fig_name, rebuild)
        # ISSUE: should consolidate this with layer.scale attribute
        # and move all to build_layer method?
        for figName, fig in self.figs.items():
            f = plt.figure(fig.fignum)
            for pos in fig.arrange.keys():
                ax = fig.arrange[pos]['axes_obj']
                xdom, ydom = fig.arrange[pos]['scale']
                ax.set_xlim(xdom)
                ax.set_ylim(ydom)
            f.canvas.draw()
        if not self.shown:
            # TEMP
            plt.ion() #Artists not appearing on axes without call to ion
            #plt.show()
            self.shown = True

        do_save = False # default

        if ignore_wait:
            wait = False
        else:
            wait = self.wait_status
        # force wait overrides, if set
        if force_wait is not None:
            wait = force_wait
        if wait:
            print("Commands:\n=========")
            print("N <RETURN>: Stop waiting on each iteration")
            print("A <RETURN>: Stop waiting and save all figures on iterations")
            print("S <RETURN>: Save this figure and continue")

            key = input('Enter command or <RETURN> to continue or ^D to quit: ')
            if key in ['N', 'n']:
                self.wait_status = False
            elif key in ['S', 's']:
                do_save = True
            elif key in ['A', 'a']:
                self.save_status = True
                self.wait_status = False

        if self.save_status or do_save:
            self.gui.save(None)
            #if self.dm is not None:
                #dirpath = self.dm._dirpath
            #else:
                #dirpath = ''
            #f.savefig(os.path.join(dirpath, get_unique_name(fig_name,
                                                            #start=1)+'.png'),
                      #format='png')


    def build_layers(self, layer_list, ax, rescale=None, figure=None, rebuild=False):
        """
        Convenience function to group layer refresh/build calls. Current figure for
        given list of layer names is assumed unless optionally specified.

        Optional rebuild = True argument will rebuild all plots, overwriting any
        previous object handles.
        """
        fig_struct, figure = self._resolve_fig(figure)
        if not fig_struct.display:
            return

        for layer_name in layer_list:
            self.build_layer(figure, layer_name, ax, rescale, force=(rebuild or
                                                                    fig_struct['layers'][layer_name].force))


    def update_dynamic(self, time, dynamicFns, hard_reset=False):
        """
        Dynamic callback functions always accept time as first argument.
        Optional second argument is hard_reset Boolean.
        """
        for figName in self.figs:
            for layer in self.figs[figName].layers:
                argsLayer = self.figs[figName].layers[layer]

                if not argsLayer.dynamic:
                    continue
                else:
                    assert isinstance(dynamicFns, dict), \
                           "Dynamic functions must be a dictionary of layer-function pairs"

                    if layer in dynamicFns:
                        #print("update_dynamic calling function: %s" % str(dynamicFns[layer]))
                        dynamicFns[layer](time, hard_reset)


    def _resolve_fig(self, figure):
        """
        Internal utility to return a figure structure and figure name,
        given that the figure argument may be None (selecting the current figure)
        """
        if figure is None and self.currFig is None:
            raise ValueError("Must set current figure")
        elif figure is None and self.currFig is not None:
            figure = self.currFig

        try:
            fig_struct = self.figs[figure]
        except KeyError:
            raise ValueError("Figure %s does not exist"%figure)

        return fig_struct, figure


    def _resolve_layer(self, figure, layer):
        """
        Internal utility to find the named layer in the given figure, and return its
        struct.
        """
        try:
            fig_struct = self.figs[figure]
        except KeyError:
            raise KeyError("Invalid figure name: %s" % figure)
        try:
            return fig_struct.layers[layer]
        except KeyError:
            raise KeyError("Invalid layer name: %s in figure %s" % (layer, figure))

    def _retrive_subplots(self, layer):
        """
        Internal utility to find all subplots a given layer has been assigned to through arrange_fig.
        """
        subplots = []
        for sp, dic in self.figs[self.currFig]['arrange'].items():
            if dic['layers'] is '*':
                subplots = list(self.figs[self.currFig]['arrange'].keys())
                break

            if layer in dic['layers']:
                subplots += [sp]

        return subplots


    def build_layer(self, figure_name, layer_name, ax, rescale=None, force=False):
        """
        Consolidates layer information into matplotlib.artist objects
        rescale (pair of pairs) may be set if the axes' current scalings
           should be preserved, overriding the layer's set scale, if any

        Optional force = True argument will rebuild plot object handles.
        """
        fig_struct = self.figs[figure_name]
        if not fig_struct.display:
            return

        lay = fig_struct.layers[layer_name]
        if not lay.display:
            # switch off visibility of all object handles
            for obj in lay.handles.values():
                obj.set_visible(False)
            return

        ##ISSUE: This is very finnicky. Should find more reliable way to clear artists in matplotlib axes.
        if force:
            for h in lay.handles.values():
                ##ISSUE: Error in MPL: list.remove(x): x not in list, when using HH_simple_demo
                try:
                    #ax.lines.remove(h)
                    h.remove()
                except ValueError:
                    pass

            lay.handles = OrderedDict({})

        ##ISSUE: Would be preferable to only do this if we know these fields have been changed.
        for hname, handle in lay.handles.items():
            ##ISSUE: Sometimes an hname isn't in the lay.data. Should not have to use a try except here.
            try:
                handle.set_linewidth(lay.data[hname]['linewidth'])
                handle.set_markersize(lay.data[hname]['markersize'])
                handle.set_zorder(lay.data[hname]['zorder'])
                handle.set_color(lay.data[hname]['style'][0])
            except (KeyError, AttributeError) as e:
                pass

        for dname, dstruct in lay.data.items():
            # we should have a way to know if the layer contains points
            # that may or may not already be updated *in place* and therefore
            # don't need to be redrawn

            ##ISSUE: This block doesn't seem to do anything.
            try:
                obj = lay.handles[dname]
            except KeyError:
                # object not actually plotted yet
                pass
            else:
                obj.set_visible(dstruct['display'])

            # For now, default to first subplot with 0 indexing if multiple exist
            if dstruct['subplot'] == None:
                dstruct['subplot'] = self._retrive_subplots(layer_name)[0]

            #Use subplot string for current data to retrieve the axes object.
            ax = self.figs[self.currFig].arrange[dstruct['subplot']]['axes_obj']

            try:
                # process user-defined style
                s = dstruct['style']

                if isinstance(s, str):
                    style_as_string = True
                elif isinstance(s, dict):
                    style_as_string = False
                else:
                    style_as_string = False

                if s == "" or s is None:
                    # default to black lines
                    s = 'k-'
            except KeyError:
                pass

            # in case in future we wish to have option to reverse axes
            ix0, ix1, ix2 = 0, 1, 2

            #if dstruct['selected']:
                #linewidth = 2.5
                #markersize = 12
            #else:
                #linewidth = 1
                #markersize = 6

            if lay.kind == 'text' and dstruct['display']:
                ##ISSUE: These probably just need to be "if dname not in handles:", since force = True
                ## will clear lay.handles anyways.
                if dname not in lay.handles or force:
                    if dstruct['use_axis_coords']:
                        lay.handles[dname] = ax.text(dstruct['data'][ix0],
                             dstruct['data'][ix1],
                             dstruct['text'],
                             transform=ax.transAxes,
                             fontsize=20, color=s[0])
                    else:
                        lay.handles[dname] = ax.text(dstruct['data'][ix0],
                             dstruct['data'][ix1],
                             dstruct['text'],
                             fontsize=20, color=s[0])

            elif lay.kind == 'patch':
                pos = dstruct['data']
                for i in range(len(pos[0])):
                    #This must generalize to other patches.
                    lay.handles[dname] = ax.add_artist(dstruct['patch']((pos[0][i], pos[1][i]),
                                  dstruct['radius'][i],
                                  color = dstruct['color'],
                                  visible = dstruct['display']))

            elif lay.kind == 'obj':

                coords = dstruct['data']
                try: #Line
                    l = dstruct['obj'](coords[0], coords[1], linewidth=dstruct['linewidth'], color='y', visible= dstruct['display'])
                except TypeError: #Rectangle
                    l = dstruct['obj'](coords[0], coords[1][0], coords[1][1], linewidth=dstruct['linewidth'], color='y', visible= dstruct['display'], fill= False)

                ##ISSUE: try/except not needed here.
                try:
                    self.gui.context_objects[dname].handle = l
                    lay.handles[dname] = ax.add_artist(l)
                    lay.handles[dname].set_picker(2.5)
                except KeyError:
                    pass

            elif lay.kind == 'data':
                if dname not in lay.handles or force:
                    try:
                        lay.handles[dname] = ax.add_collection(dstruct['data'])
                    except AttributeError:
                        if style_as_string:
                            #Check if data are two or three dimensional.
                            if len(dstruct['data']) == 2:
                                ##ISSUE: Should repeat changes made to this case (i.e., setting the picker and new dstruct properties)
                                ## to the other cases. At least the other 2d case.
                                ## style_as_string == False is a potential landmine.
                                lay.handles[dname] = \
                                    ax.plot(dstruct['data'][ix0], dstruct['data'][ix1],
                                            s, linewidth= dstruct['linewidth'],
                                            zorder = dstruct['zorder'], markersize = dstruct['markersize'],
                                            visible= dstruct['display'])[0]
                                #ax.add_artist(lay.handles[dname])
                                lay.handles[dname].set_picker(2.5)

                            elif len(dstruct['data']) == 3:
                                lay.handles[dname] = \
                                    ax.plot(dstruct['data'][ix0], dstruct['data'][ix1], dstruct['data'][ix2],
                                            s, visible= dstruct['display'])[0]
                        else:
                            #Display? Visibility?
                            if len(dstruct['data']) == 2:
                                lay.handles[dname] = \
                                    ax.plot(dstruct['data'][ix0], dstruct['data'][ix1],
                                            **s)[0]
                                #ax.add_artist(lay.handles[dname])
                                lay.handles[dname].set_picker(True)
                            elif len(dstruct['data']) == 3:
                                lay.handles[dname] = \
                                    ax.plot(dstruct['data'][ix0], dstruct['data'][ix1], dstruct['data'][ix2],
                                            **s)[0]

                    #ax.add_artist(lay.handles[dname])
                    lay.force = False

        if rescale is not None:
            # overrides layer scale
            sc = rescale
        elif lay.scale is not None:
            sc = lay.scale
        else:
            sc = None
        if sc is not None:
            try:
                ax.set_xlim(sc[0])
            except TypeError:
                pass
            try:
                ax.set_ylim(sc[1])
            except TypeError:
                pass


class diagnosticGUI(object):
    """
    Diagnostics Graphical User Interface
    """
    def __init__(self, objPlotter, points=None, verbose_level=1):

        assert isinstance(objPlotter, Plotter), \
               "Plotter_GUI must be instantiated with a Plotter object."

        self.cb_axes = []

        # easy user access to these basic attributes: current time and
        # index into self.points
        self.t = None
        self.ix = None

        self.points = None
        # ISSUE: add Trajectory object too?

        # times is an array from trajectory points
        self.times = None

        # clipboard for a point selected in a dynamic sub-plot (e.g., phaseplane)
        self.clipboardPt = None

        # Captured points in all sub-plots at current time, from GUI button
        self.capturedPts = {}

        # internal verbosity level
        self.verbose = verbose_level

        if points is not None:
            self.add_points(points)

        # ---------------
        # Rocket stuff
        # ---------------

        self.model = None
        self.gen_versioner = None #Universal
        self.context_changed = None #Universal


        # ---------------
        # Internal stuff
        # ---------------

        # masterWin is the figure handle for main GUI window
        self.masterWin = None

        # default: does not expect time-parameterized trajectories
        # in main window
        self._with_times = False

        self.timePlots = []
        self.timeLines = []

        self.plotter = objPlotter
        self.plotter.gui = self

        self.widgets = {}

        # callback functions for dynamic plots indexed by layer name
        self.dynamicPlotFns = {}

        # axes objects for dynamic plots indexed by layer name
        self.dynamicPlots = {}

        self._mouseUp = True
        self._mouseDrag = False
        self._key_mod = None
        # key modifiers for +/- dt control
        self._key_mod_dix = {None: 5,
                             'shift': 10,
                             'control': 1}
        self._last_ix = None

        #Delete these.

        #global plotter
        #plotter = objPlotter

        #global gui
        #gui = self

    def initialize_callbacks(self, fig):
        #INIT FROM GUIROCKET
        self.fig = fig
        self.context_objects = {}

        self.selected_object = None
        self.selected_object_temphandle = None
        self.current_domain_handler = dom.GUI_domain_handler(self)
        self.tracked_objects = []
        self.calc_context = None

        self.last_output = None

        self.mouse_wait_state_owner = None

        ##Handled in plotter2d._subplots
        #self.RS_line = RectangleSelector(self.ax, self.onselect_line, drawtype='line')
        #self.RS_line.set_active(False)
        self.RS_boxes = {}
        self.RS_lines = {}

        self.fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
        self.fig.canvas.mpl_connect('pick_event', self.pick_on)

        evKeyOn = self.fig.canvas.mpl_connect('key_press_event', self.key_on)
        evKeyOff = self.fig.canvas.mpl_connect('key_release_event', self.key_off)

    def add_data_traj(self, traj, points=None):
        """
        Provide the trajectory (or other curve object) for the
        data to be investigated. In case trajectory is not defined by
        a standard mesh, a user-defined sampling of points can be
        optionally provided.
        """
        # ISSUE: This structure assumes time-dependent data only
        #  (not sufficiently general purpose)
        self.traj = traj
        if points is None:
            self.points = traj.sample()
        else:
            self.points = points
        try:
            self.times = self.points['t']
        except KeyError:
            # trajectory is not parameterized by 't'
            self.times = None

    def add_data_points(self, data, figure=None, layer=None, subplot=None,
                           style=None, linewidth = 1, name=None, display=True,
                           force=False, log=None, coorddict=None):
        maxspeed = 2.2 #ISSUE: This should be replaced with something general purpose. Borrowed from Bombardier.

        try:
            fig_struct, figure = self.plotter._resolve_fig(None)
        except ValueError:
            self.plotter.clean()
            self.plotter.add_fig('Master', domain= [(0,1), (0,1)])
            fig_struct, figure = self.plotter._resolve_fig(None)

        if isinstance(data, numpy.ndarray) or isinstance(data, list):
            self.plotter.add_data(data, figure=figure, layer=layer, subplot=subplot,
                           style=style, name=name,
                           display=display,
                           linewidth=linewidth,
                           force=force, log=log)

        elif isinstance(data, Points.Pointset) or isinstance(data, pp.Point2D):
            if coorddict is None:
                warnings.warn("coorddict needed to specify how Pointset coords should be plotted.")
                return

            else:
                addingDict = {}

                for key, val in coorddict.items():
                    #Extract x and y data.
                    try:
                        xs = data[val.get('x')]
                    except (IndexError, StopIteration) as e:
                        try:
                            del xs
                        except UnboundLocalError:
                            pass

                    try:
                        ys = data[val.get('y')]
                    except (IndexError, StopIteration) as e:
                        try:
                            del ys
                        except UnboundLocalError:
                            pass

                    if not key in addingDict.keys() and ('xs' in locals() or 'ys' in locals()):
                        addingDict[key] = {}

                    try:
                        addingDict[key]['data'] = [xs, ys]

                    #If only x/y provided, use key data for other coordinate.
                    except NameError:
                        try:
                            addingDict[key]['data'] = [xs, data[key]]
                        except UnboundLocalError:
                            pass
                        try:
                            addingDict[key]['data'] = [data[key], ys]
                        except UnboundLocalError as e:
                            pass

                    #Extract object
                    if val.get('object') == 'collection':
                        addingDict[key]['segments'] = [( (xs[i], ys[i]), (xs[i+1], ys[i+1]) ) for i in range(len(xs)-1)]
                    elif val.get('object') == 'circle':
                        addingDict[key]['patch'] = plt.Circle

                    #Extract style
                    try:
                        addingDict[key]['style'] = val['style']
                    except KeyError:
                        pass

                    #Extract layer
                    try:
                        if val['layer'] not in fig_struct.layers.keys():
                            self.plotter.add_layer(val['layer'])
                            print("Added new layer %s to plotter."%val['layer'])
                        addingDict[key]['layer'] = val['layer']
                    except KeyError:
                        pass

                    #Extract name
                    try:
                        addingDict[key]['name'] = val['name']
                    except KeyError:
                        pass

                    #Extract radii
                    try:
                        addingDict[val['map_radius_to']]['radius'] = data[key]
                    except:
                        try:
                            addingDict[val['map_radius_to']] = {}
                            addingDict[val['map_radius_to']]['radius'] = data[key]
                        except KeyError:
                            pass

                    #Perform color mapping
                    try:
                        vals = data[key]
                        norm = mpl.colors.Normalize(vmin=0, vmax=maxspeed)
                        cmap=plt.cm.jet #gist_heat
                        try:
                            addingDict[val['map_color_to']]['style'] = cmap(norm(vals))
                        except KeyError:
                            addingDict[val['map_color_to']] = {}
                            addingDict[val['map_color_to']]['style'] = cmap(norm(vals))

                    except KeyError:
                        pass
                    except IndexError:
                        pass

                #add_data for each plotting variable in the pointset.
                for key, val in addingDict.items():
                    try:
                        linecollection = mpl.collections.LineCollection(val['segments'], colors=val['style'])
                        #addingDict[key]['traj'] = self.reducePointset(data, coorddict, list(addingDict.keys())[0])
                        addingDict[key]['data'] = linecollection
                    except KeyError:
                        #addingDict[key]['traj'] = self.reducePointset(data, coorddict, list(addingDict.keys())[0])
                        pass

                    try:
                        tra = self._reduce_pointset(data, coorddict, list(addingDict.keys())[0])
                    except (KeyError, ValueError) as e:
                        tra = None

                    try:
                        lay = addingDict[key]['layer']
                    except KeyError:
                        lay = None

                    try:
                        nam = addingDict[key]['name']
                    except KeyError:
                        nam = None

                    try:
                        self.plotter.add_patch(addingDict[key]['data'], addingDict[key]['patch'],
                                         layer = lay,
                                         name = nam,
                                         force = True,
                                         radius = addingDict[key]['radius'],
                                         color = addingDict[key]['style'])

                    except:
                        self.plotter.add_data(addingDict[key]['data'],
                                        style = addingDict[key]['style'],
                                        layer = lay,
                                        name = nam,
                                        traj = tra,
                                        linewidth = linewidth, ##ISSUE: Should do this through coorddict as well.
                                        force = True)
        else:
            warnings.warn("add_data_points received an unsupported type for parameter data")
            return

    def _reduce_pointset(self, ptset, coorddict, key):
        """
        Convenience function for extracting two desired points from a pointset
        with many variables.
        """
        try:
            x = coorddict[key]['x']
        except KeyError:
            x = key

        try:
            y = coorddict[key]['y']
        except KeyError:
            y = key

        try:
            new_pts = Pointset({'coordarray': [ptset[x], ptset[y]],
                         'coordnames': ['x', 'y'],
                         'indepvarname':'t',
                         'indepvararray': ptset.indepvararray})
        except (AttributeError, TypeError) as e:
            new_pts = Pointset({'coordarray': [ptset[x], ptset[y]],
                                'coordnames': ['x', 'y']})
        return new_pts


    def add_widget(self, widg, axlims, callback=None, **kwargs):
        """
        Create a matplotlib widget
        """
        if not issubclass(widg, mpl.widgets.Widget):
            raise TypeError("widg must be a subclass of matplotlib.widget")
        if not callable(callback) and callback is not None:
            raise TypeError("callback must be callable")

        ax=plt.axes(axlims)

        widget = widg(ax=ax, **kwargs)
        self.widgets[kwargs['label']] = widget
        try:
            self.widgets[kwargs['label']].on_changed(callback)
        except AttributeError:
            self.widgets[kwargs['label']].on_clicked(callback)


    def add_time_from_points(self, points):
        self.traj = None
        self.points = points
        try:
            self.times = points['t']
        except KeyError:
            self.times = None
        except PyDSTool_KeyError:
            pass


    def build_plotter(self, figsize=None, with_times=True, basic_widgets=True, callbacks_on= True):
        """
        Create time bar widget.
        Create capture points widget.
        Also closes all figure windows and resets masterWin attribute.

        Optional size is a pair of figure screen size measurements in inches
        """

        plt.close('all')
        self.masterWin = None
        self._with_times = with_times

        for figName, fig_struct in self.plotter.figs.items():
            if figsize is not None:
                fig_handle = plt.figure(fig_struct.fignum, figsize=figsize)
                # fig sizing later, with fig.set_figsize_inches doesn't seem
                # to work properly unless put in the initialization call
            else:
                fig_handle = plt.figure(fig_struct.fignum)
            fig_handle.canvas.set_window_title(fig_struct.title + " : Master window")
            if figName != 'Master' and figName != 'master':
                continue
            # ====== Set up master window controls
            plt.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.1,
                               wspace=0.2, hspace=0.23)

            self.masterWin = fig_handle

            if callbacks_on:
                self.initialize_callbacks(self.masterWin)

            # Time bar controls time lines in figures
            # ISSUE: Not all uses of this class use time
            if with_times:
                sliderRange = self.times
                slide = plt.axes([0.25, 0.02, 0.65, 0.03])
                tMin = min(sliderRange)
                tMax = max(sliderRange)
                if self.t is None:
                    self.set_time( (tMin + tMax)/2. )
                self.widgets['timeBar'] = Slider(slide, 'Time', tMin, tMax,
                                                valinit=self.t, color='r',
                                                dragging=False, valfmt='%1.4f')

                # button axes are in figure coords: (left, bottom, width, height)

                ## +/- dt buttons
                m_dt_Button = Button(plt.axes([0.16, 0.02, 0.017, 0.03]), '-dt')
                self.widgets['minus_dt'] = m_dt_Button

                p_dt_Button = Button(plt.axes([0.18, 0.02, 0.017, 0.03]), '+dt')
                self.widgets['plus_dt'] = p_dt_Button

            if basic_widgets:
                # Capture point button in lower left
                captureButton = Button(plt.axes([0.055, 0.02, 0.08, 0.03]), 'Capture Point')
                self.widgets['capture_point'] = captureButton

                # Refresh button
                refreshButton = Button(plt.axes([0.005, 0.02, 0.045, 0.03]), 'Refresh')
                self.widgets['refresh'] = refreshButton

                # Go back to last point button
                backButton = Button(plt.axes([0.005, 0.06, 0.045, 0.03]), 'Back')
                self.widgets['go_back'] = backButton

                # Go back to last point button
                saveButton = Button(plt.axes([0.055, 0.06, 0.08, 0.03]), 'Save')
                self.widgets['save'] = saveButton

                # Display graphics objects hierarchy
                showTreeButton = Button(plt.axes([0.86, 0.02, 0.12, 0.03]), 'Show Tree')
                self.widgets['showTree'] = showTreeButton

                self.widgets['capture_point'].on_clicked(self.capture_point)
                self.widgets['refresh'].on_clicked(self.refresh)
                self.widgets['go_back'].on_clicked(self.go_back)
                self.widgets['save'].on_clicked(self.save)
                self.widgets['showTree'].on_clicked(self.show_tree)

            self.plotter.show(update='all', rebuild=True, force_wait=False)

            # Build up each subplot, left to right, top to bottom
            shape = fig_struct.shape
            for i in range(shape[0]):
                for j in range(shape[1]):
                    ixstr = str(i+1) + str(j+1)
                    try:
                        subplot_struct = fig_struct.arrange[ixstr]
                    except (KeyError, TypeError):
                        # type error if arrange is empty list (e.g. for shape=[1,1])
                        continue
                    layer_info = subplot_struct['layers']
                    if not isinstance(layer_info, list):
                        if layer_info == '*':
                            layer_info = list(fig_struct.layers.keys())
                        else:
                            # singleton string layer name
                            layer_info = [layer_info]

                    ax = subplot_struct['axes_obj']

                    for layName in layer_info:
                        if layName in self.dynamicPlotFns:
                            self.dynamicPlots[layName] = ax
                            # initialize the layer's dynamic stuff
                            self.dynamicPlotFns[layName](self.t)
                    if with_times:
                        # add vertical time line in all time plots
                        # (provide option for user to specify as time or t)
                        if subplot_struct['axes_vars'][0].lower() in ["time", 't']:
                            self.timeLines.append(ax.axvline(x=self.t, color='r',
                                                         linewidth=3, linestyle='--'))
                            self.timePlots.extend(layer_info)

        if [isinstance(fig_struct['arrange'][pos]['axes_obj'], Axes3D) for pos in fig_struct['arrange'].keys()]:
            print("3D Axes can be rotated by clicking and dragging.")

        # Activate button & slider callbacks
        #self.widgets['capture_point'].on_clicked(self.capture_point)
        #self.widgets['refresh'].on_clicked(self.refresh)
        #self.widgets['go_back'].on_clicked(self.go_back)
        #self.widgets['save'].on_clicked(self.save)
        if with_times:
            self.widgets['timeBar'].on_changed(self.update_plots)
            self.widgets['minus_dt'].on_clicked(self.minus_dt)
            self.widgets['plus_dt'].on_clicked(self.plus_dt)

        # Activate general mouse click callbacks
        evMouseDown = fig_handle.canvas.mpl_connect('button_press_event', self.mouse_down)
        evMouseUp = fig_handle.canvas.mpl_connect('button_release_event', self.mouse_up)
        evMouseMove = fig_handle.canvas.mpl_connect('motion_notify_event', self.mouse_move)
        evKeyOn = fig_handle.canvas.mpl_connect('key_press_event', self.modifier_key_on)
        evKeyOff = fig_handle.canvas.mpl_connect('key_release_event', self.modifier_key_off)

    def clear_data(self, layer, data_name=None):
        """
        Delete data in a given layer. If no name for data is given, all data in that layer are deleted.

        ISSUE: Consider merging with clear_axes.
        """
        fig_struct, figure = self.plotter._resolve_fig(None)
        layer_struct = self.plotter._resolve_layer(figure, layer)

        if data_name is not None:
            del layer_struct.data[data_name]
        else:
            layer_struct.data = {}

    def clear_axes(self, subplot, figure=None):
        """
        Clears lines and points in sub-plot axes (given as an axis
        object or sub-plot string name) without removing
        other items such as title, labels.
        """
        fig_struct, figure = self.plotter._resolve_fig(figure)
        arrPlots = fig_struct.arrange
        subplot_struct = None
        if isinstance(subplot, str):
            for ixstr, spec in arrPlots.items():
                if subplot in [ixstr, spec['name']]:
                    subplot_struct = fig_struct.arrange[ixstr]
                    ax = subplot_struct['axes_obj']
        else:
            # try subplot as axes object
            try:
                ixstr = self.plotter.subplot_lookup[subplot][2]
            except IndexError:
                pass
            else:
                ax = subplot
                subplot_struct = fig_struct.arrange[ixstr]

        if subplot_struct is None:
            raise ValueError("subplot argument not a valid index string, name, or axes object")

        # ax.clear correctly removes some of the plots that persist
        # with ax.lines = [], but it also clears the axes labels
        ax.lines = []
        #ax.clear()
##        ax.set_title(subplot_struct['name'])
##        axes_vars = subplot_struct['axes_vars']
##        ax.set_xlabel(axes_vars[0])
##        ax.set_ylabel(axes_vars[1])
##        ax.figure.canvas.draw()


    # ======== Standard callback functions for the GUI display

    def mouse_down(self, ev):
        self._mouseUp = False
        #print("mouse_down", self._mouseUp)

    def mouse_up(self, ev):
        # NOTE: zoom dragging will not get completed before this callback
        # so trying to refresh as a result of zoom here will fail
        self._mouseUp = True
        #print("mouse_up", self._mouseUp)
        # if in a time-based sub-plot and not dragging
        do_get = False
        if not self._mouseDrag:
            do_get = True and self._with_times
            # check widget axes
            for w in self.widgets.values():
                if ev.inaxes == w.ax:
                    do_get = False
        if do_get:
            # resolve whether up happens in time-based plot
            # or a dynamic plot
            if ev.inaxes in self.dynamicPlots.values():
                self.get_dynamic_point(ev)
            else:
                self.get_point(ev)
        self._mouseDrag = False


    def mouse_move(self, ev):
        if self._mouseUp:
            self._mouseDrag = False
        else:
            self._mouseDrag = True

    def add_legend(self, colors, labels, subplot):
        """
        Creates a matplotlib legend associated with a given subplot.

        May not be compatible with versions earlier than python 2.7 and matplotlib 1.2.
        """
        fig_struct, fig = self.plotter._resolve_fig(None)

        if not isinstance(colors, list) or not isinstance(labels, list):
            raise TypeError("colors and labels must be lists.")

        if len(colors) != len(labels):
            raise ValueError("colors and labels lists must be the same length.")

        handles = []
        for i in range(len(colors)):
            handles.append(mpl.patches.Patch(color= colors[i], label= labels[i]))

        fig_struct.arrange[subplot]['legend'] = handles


    def get_point(self, ev):
        """
        If mouse click is released inside a time-based plot, change current time
        on time bar to that selected.
        """
        self.set_time(ev.xdata)

    def set_point(self, name, pt, layer, figure=None):
        """
        Wrapper method for plotter.set_point
        """
        self.plotter.set_point(name, pt, layer)

    def add_layer(self, layer_name, figure=None, set_to_active=True, subplot=None, **kwargs):
        """
        Wrapper method for plotter.add_layer
        """
        self.plotter.add_layer(layer_name, figure=figure, set_to_active=set_to_active, subplot=subplot, **kwargs)

    def add_fig(self, label, title="", xlabel="", ylabel="", tdom=None, domain=None, display=True):
        """
        Wrapper method for plotter.add_fig
        """
        self.plotter.add_fig(label, title=title, xlabel=xlabel, ylabel=ylabel, tdom=tdom,
                           domain=domain, display=display)

    def show(self, update='current', rebuild=False, force_wait=None):
        """
        Wrapper method for plotter.show
        """
        self.plotter.show(update= update, rebuild= rebuild, force_wait= force_wait)

    def show_legends(self, figure=None, subplot=None):
        """
        Wrapper method for plotter.show_legends
        """
        self.plotter.show_legends(figure=figure, subplot=subplot)


    def clean(self):
        """
        Wrapper method for plotter.clean
        """
        self.plotter.clean()

    def build_layers(self, layer_list, ax, rescale=None, figure=None,
                    rebuild=False):
        """
        Wrapper method for plotter.build_layers
        """
        self.plotter.build_layers(layer_list, ax, rescale=rescale, figure=figure,
                                rebuild=rebuild)

    def get_dynamic_point(self, ev):
        """
        If mouse clicks inside a user-specified 'dynamic' sub-plot axis,
        put the (x,y) coords of the click as a point onto the clipboard.
        """
        figName, layer_info, ixstr = self.plotter.subplot_lookup[ev.inaxes]
        fig_struct = self.plotter.figs[figName]
        # all layers share the same axes, so just get the first
        layer_struct = fig_struct.layers[layer_info[0]]
        self.clipboardPt = pp.Point2D(ev.xdata, ev.ydata,
                                   xname=layer_struct.axes_vars[0],
                                   yname=layer_struct.axes_vars[1])
        print("Clipboard now contains: %s" % str(self.clipboardPt))


    def capture_point(self, ev):
        """
        Create a dict of Point objects from the current time point in all sub-plots
        of all figures. Keys of the dict are the figure names
        Stored in the attribute capturedPts

        For non time-based data, capture_point creates a dictionary of related data_GUI objects.
        Data are considered related if they share the same name (despite being in different layers).
        """
        pts_dict = {}
        # pts_dict : figName -> Point with coords t and <layer_names>
        # layer_dict : layer_name -> value

        for figName in self.plotter.figs:

            fig_struct = self.plotter.figs[figName]


            if self.t is not None:
                pt_dict = {'t': self.t, 'ix': self.ix}

                for subplot, subplot_struct in fig_struct.arrange.items():
                    layer_info = subplot_struct['layers']

                    for layName in layer_info:
                        #yvar = fig_struct.layers[layName].axes_vars[1]
                        # ignore dynamic sub-plots such as phase plane, which don't show
                        # time-varying quantities that can be sampled this way
                        if layName in self.timePlots:
                            for data_name, traj in fig_struct.layers[layName].trajs.items():
                                if fig_struct.layers[layName].kind == 'data':
                                    pt_dict[data_name] = traj(self.t)['y'] #ISSUE: Why is it hardwired to select y coord?

            elif self.selected_object is not None:
                pt_dict = {}

                for layer, layer_struct in fig_struct.layers.items():
                    for dname, dstruct in layer_struct.data.items():

                        if dname == self.selected_object.name:
                            key = dname+'_'+layer
                            pt_dict[key] = data_GUI(dname, layer_struct.handles[dname], layer, self)

            else:
                print("No selected object or points at time-step to capture")
                return

            print("figName: ")
            print(figName)
            print("pt_dict: ")
            print(pt_dict)
            print(type(pt_dict))

            if self.t is not None:
                print("Point(pt_dict): ")
                print(Point(pt_dict))
                type(Point(pt_dict))
                print("...")
                pts_dict.update({figName: Point(pt_dict)})

            else:
                print("...")
                pts_dict.update({figName: pt_dict})

        if self.verbose >= 1 and self.t is not None:
            # print point information to stdout console
            for figName, pt in pts_dict.items():
                print("\n\n==============================")
                print("Figure: %s" % figName)
                print("@ time = %.3f" % pt['t'])
                print(pt)
            print("\n")

        self.capturedPts = pts_dict


    def set_time(self, time):
        """
        Find nearest time in data set to selected value, and set attribute t
        accordingly. Also updates time bar widget if it has been
        created.
        """
        if time is None:
            # do nothing at all on mouse drag event in window
            return
        ix = int(np.argmin(abs(self.times - time)))
        if ix == self.ix:
            # nothing to do
            return
        self.t = self.times[ix]
        self._last_ix = self.ix
        self.ix = ix
        do_draw = False
        for line in self.timeLines:
            line.set_data(([time, time], [0,1]))
            do_draw = True
        try:
            if self.widgets['timeBar'].val != time:
                self.widgets['timeBar'].set_val(time)
        except KeyError:
            # timeBar not yet created
            pass
        else:
            do_draw = True
        if do_draw:
            plt.draw()

    def go_back(self, ev):
        if self._last_ix is not None:
            self.update_plots(self.times[self._last_ix])

    def update_plots(self, new_time):
        self.set_time(new_time)
        self.plotter.update_dynamic(self.t, self.dynamicPlotFns)

    def refresh(self, ev):
        """
        For refresh button, e.g. use after zoom in dynamic plot
        """
        hard_reset = self._key_mod == 'shift'
        self.plotter.update_dynamic(self.t, self.dynamicPlotFns,
                                   hard_reset)

    def save(self, ev):
        """
        For save button. Saves current figure as a .png
        in working directory or dm directory if provided.
        """
        fig_struct, fig_name = self.plotter._resolve_fig(self.plotter.currFig)
        f = plt.figure(fig_struct.fignum)

        if self.plotter.dm is not None:
            dirpath = self.plotter.dm._dirpath
        else:
            dirpath = ''
        f.savefig(os.path.join(dirpath, get_unique_name(fig_name,
                                                        start=1)+'.png'),
                  format='png')

    def minus_dt(self, ev):
        """For the -dt button.
        Moves -1 index position; holding shift moves -10 indices
        """
        dix = self._key_mod_dix[self._key_mod]
        ix = max(0, min(self.ix-dix, len(self.times)-1))
        if ix != self.ix:
            self.set_time(self.times[ix])
        # else do nothing

    def plus_dt(self, ev):
        """For the +dt button.
        Moves +1 index position; holding shift moves +10 indices
        """
        dix = self._key_mod_dix[self._key_mod]
        ix = max(0, min(self.ix+dix, len(self.times)-1))
        if ix != self.ix:
            self.set_time(self.times[ix])
        # else do nothing

    def modifier_key_on(self, ev):
        self._key_mod = ev.key

    def modifier_key_off(self, ev):
        self._key_mod = None

    def pick_on(self, event):
        """
        Pick artists in axes and set them as the selected object.
        """
        ##ISSUE: If click is within range of multiple artist, pick_on is called many times in succession,
        ##and picks all those artists as selected objects in sequence.
        fig_struct, fig = self.plotter._resolve_fig(None)
        is_con_obj = False

        if isinstance(event.artist, mpl.lines.Line2D):
            artist_data = event.artist.get_data()

        elif isinstance(event.artist, mpl.patches.Rectangle):
            #Need to use width of rectangle to CALCUlATE other vertices.
            artist_data = [[event.artist.get_x(), event.artist.get_x()+ event.artist.get_width()],
                           [event.artist.get_y(), event.artist.get_y()+ event.artist.get_height()]]

        for con_obj in self.context_objects.values():
            if isinstance(con_obj, shape_GUI) and \
               artist_data[0][0] == con_obj.x1 and \
               artist_data[0][1] == con_obj.x2 and \
               artist_data[1][0] == con_obj.y1 and \
               artist_data[1][1] == con_obj.y2:
                self.set_selected_object(con_obj)

                is_con_obj = True

            if isinstance(con_obj, domain_GUI):
                ##ISSUE: domain_GUI picking not yet implemented.
                is_con_obj = True
                pass


        if not is_con_obj:
            for subplot, subplot_struct in fig_struct.arrange.items():
                if event.mouseevent.inaxes is subplot_struct['axes_obj']:
                    for lay in subplot_struct['layers']:
                        layer_struct = self.plotter._resolve_layer(fig, lay)
                        for name, artist in layer_struct.handles.items():
                            if artist is event.artist:
                                self.set_selected_object(data_GUI(name, artist, lay, self))

        self.user_pick_func(event)


    def show_tree(self, event= None):
        """
        Prints to terminal the structure of the current Fovea figures and all the graphical objects
        (layers, context objects and data) contained therein.
        """
        for fig_name, fig_struct in self.plotter.figs.items():
            print('Figure: %s'%fig_name)
            for lay_name, lay_struct in fig_struct['layers'].items():
                print('--Layer: %s'%lay_name)
                for data_name, data_struct in lay_struct['data'].items():
                    print('----%s: %s' %(lay_struct.kind, data_name))

    def navigate_selected_object(self, k):
        """
        Key presses used to manipulate the currently selected object.
        """
        so = self.selected_object
        fig_struct, figure = self.plotter._resolve_fig(None)
        layer_struct = self.plotter._resolve_layer(self.plotter.currFig, so.layer)

        #If context object.
        if isinstance(so, line_GUI) or isinstance(so, box_GUI):

            ax = fig_struct['arrange'][layer_struct['data'][so.name]['subplot']]['axes_obj']

            xl = ax.get_xlim()
            yl = ax.get_ylim()

            step_sizeH = 0.02*abs(xl[0]-xl[1])
            step_sizeV = 0.02*abs(yl[0]-yl[1])
            nav = False

            if k == 'n':
                ##ISSUE: Currently causes "RuntimeError: can't re-enter readline"
                name = input('Enter new name for the selected object: ')
                so.update(name = name)

            if k == 'backspace':
                so.remove(draw= True)
                self.selected_object = None

            if isinstance(so, line_GUI):
                if k == 'left':
                    so.update(x1 = (so.x1 - step_sizeH), x2 = (so.x2 - step_sizeH))
                    #so.update(x = [so.x1 - step_sizeH, so.x2 - step_sizeH])
                    nav = True
                elif k == 'right':
                    so.update(x1 = (so.x1 + step_sizeH), x2 = (so.x2 + step_sizeH))
                    #so.update(x = [so.x1 + step_sizeH, so.x2 + step_sizeH])
                    nav = True
                elif k == 'up':
                    so.update(y1 = (so.y1 + step_sizeV), y2 = (so.y2 + step_sizeV))
                    #so.update(y = [so.y1 - step_sizeV, so.y2 - step_sizeV])
                    nav = True
                elif k == 'down':
                    so.update(y1 = (so.y1 - step_sizeV), y2 = (so.y2 - step_sizeV))
                    #so.update(y = [so.y1 + step_sizeV, so.y2 + step_sizeV])
                    nav = True

            if isinstance(so, box_GUI):
            ##ISSUE: After a box_GUI has been moved, it seems to become un-pickable.
                if k == 'left':
                    so.update(x1 = (so.x1 - step_sizeH))
                    #nav = True
                elif k == 'right':
                    so.update(x1 = (so.x1 + step_sizeH))
                    #nav = True
                elif k == 'up':
                    so.update(y1 = (so.y1 + step_sizeV))
                    #nav = True
                elif k == 'down':
                    so.update(y1 = (so.y1 - step_sizeV))
                    #nav = True

            #if nav:
                #self.user_nav_func()
                #nav = False

            if k == 'm':
                so.update(x1 = xl[0], x2 = xl[1], y1 = np.mean([so.y1, so.y2]), y2 = np.mean([so.y1, so.y2]))

                ##Extend horizontally or vertically, depending on angle of line.
                #if not -45 < so.ang_deg < 45:
                    #so.update(y1 = yl[0], y2 = yl[1], x1 = np.mean([so.x1, so.x2]), x2 = np.mean([so.x1, so.x2]))
                #else:
                    #so.update(x1 = xl[0], x2 = xl[1], y1 = np.mean([so.y1, so.y2]), y2 = np.mean([so.y1, so.y2]))

        elif isinstance(so, data_GUI):
            #Retrieve the layer_struct holding this handle.
            #for subplot, subplot_struct in fig_struct.arrange.items():
                #for layer in subplot_struct['layers']:
                    #ls = self.plotter._resolve_layer(figure, layer)
                    #for hname, handle in ls.handles.items():
                        #if self.selected_object == handle:
                            #layer_struct = ls
                            #break

            if k == 'up' or k == 'down':
                handles = list(layer_struct.handles.values())
                hnames = list(layer_struct.handles.keys())
                if k == 'up':
                    idx = (handles.index(so.handle)+1)%len(handles)
                elif k == 'down':
                    idx = (handles.index(so.handle)-1)%len(handles)
                #selected_handle = handles[idx]
                self.set_selected_object(data_GUI(hnames[idx], handles[idx], so.layer, self))
                self.user_pick_func(None)

                #for hname, handle in layer_struct.handles.items():
                    #if handle == so:
                        #handles= list(layer_struct.handles.values())
                        ##selected_handle = handles[(handles.index(handle)+1)%len(handles)]
                        #if k == 'up':
                            #selected_handle = handles[(handles.index(handle)+1)%len(handles)]
                        #elif k == 'down':
                            #selected_handle = handles[(handles.index(handle)-1)%len(handles)]
                        #self.set_selected_object(selected_handle)
                        #break

    def key_on(self, ev):
        self._key = k = ev.key  # keep record of last keypress
        # TEMP
        dom_key = '.'
        change_mouse_state_keys = ['l', 's', ' '] + [dom_key]

        if self.selected_object is not None \
           and not isinstance(self.selected_object, pp.Point2D):
            self.navigate_selected_object(k)

        #Toggle tools keys
        if self.mouse_wait_state_owner == 'domain' and \
           k in change_mouse_state_keys:
            # reset state of domain handler first
            self.current_domain_handler.event('clear')
        elif k == 'l':
            try:
                self.RS_lines[ev.inaxes].set_active(True)
            except KeyError:
                return
            print("Make a line of interest")
            self.mouse_wait_state_owner = 'line'
        elif k == 'b':
            try:
                self.RS_boxes[ev.inaxes].set_active(True)
            except KeyError:
                return
            print("Make a box of interest")
            self.mouse_wait_state_owner = 'box'
        elif k == ' ':
            print("Output of user function at clicked mouse point")
            self.mouse_cid = self.fig.canvas.mpl_connect('button_release_event', self.mouse_event_user_function)
            self.mouse_wait_state_owner = 'user_func'
        elif k == 's':
            print("Snap clicked mouse point to closest point on trajectory")
            self.mouse_cid = self.fig.canvas.mpl_connect('button_release_event', self.mouse_event_snap)
            self.mouse_wait_state_owner = 'snap'
        elif k == dom_key:
            if self.selected_object_temphandle is not None:
                self.current_domain_handler.event('reset')
            print("Click on domain seed point then initial radius point")
            # grow domain
            if self.current_domain_handler.func is None:
                print("Assign a domain criterion function first!")
                return
            else:
                # this call may have side-effects
                self.current_domain_handler.event('key')
                self.mouse_wait_state_owner = 'domain'

    def key_off(self, ev):
        # TEMP
        #pass
        self._key = None

    def onselect_line(self, eclick, erelease):
        #if eclick.inaxes not in self.cb_axes:
            #return

        if eclick.button == 1:
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata

            self.set_selected_object(line_GUI(self, pp.Point2D(x1, y1), pp.Point2D(x2, y2), subplot= eclick.inaxes), figure= self.plotter.currFig)
            print("Created line as new selected object, now give it a name")
            print("  by calling this object's .update() method with the name param")
            self.RS_lines[eclick.inaxes].set_active(False)

            self.mouse_wait_state_owner = None

    def onselect_box(self, eclick, erelease):
        #if eclick.inaxes not in self.cb_axes:
            #return

        if eclick.button == 1:
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata

            self.set_selected_object(box_GUI(self, pp.Point2D(x1, y1), pp.Point2D(x2, y2), subplot= eclick.inaxes), figure= self.plotter.currFig)
            print("Created box as new selected object, now give it a name")
            print("  by calling this object's .update() method with the name param")
            self.RS_boxes[eclick.inaxes].set_active(False)

            self.mouse_wait_state_owner = None

    def mouse_event_snap(self, ev):
        if ev.inaxes not in self.cb_axes:
            return

        fig_struct, figs = self.plotter._resolve_fig(None)
        trajs = []
        for layer_name in fig_struct['layers'].keys():
            try:
                #trajs.append(list(fig_struct['layers'][layer_name]['trajs'].values())[0])
                trajs += list(fig_struct['layers'][layer_name]['trajs'].values())
            except KeyError:
                pass

        if trajs == []:
            print("No trajectory defined")
            return
        print("\nClick: (%.4f, %.4f)" %(ev.xdata, ev.ydata))
        # have to guess phase, use widest tolerance

        found_pts = []
        #print('trajs', trajs)
        for traj in trajs:
            xname = traj.coordnames[0]
            yname = traj.coordnames[1]

            eps = 0.1
            eps = 200
            try:
                found_pt = pp.find_pt_nophase_2D(traj.sample(), pp.Point2D(ev.xdata, ev.ydata, xname= xname, yname= yname), eps=eps)
                found_pts.append(found_pt)
            except ValueError:
                pass

        if found_pts == []:
            print("No nearby point found. Try again")
            self.fig.canvas.mpl_disconnect(self.mouse_cid)
            return

        #Pick closest point from all trajectories searched.
        dist = np.inf
        for pt in found_pts:
            if math.hypot(ev.xdata - pt[2]['x'], ev.ydata - pt[2]['y']) < dist:
                data = pt
                dist = math.hypot(ev.xdata - pt[2]['x'], ev.ydata - pt[2]['y'])

        self.last_output = data
        x_snap = data[2]['x']
        y_snap = data[2]['y']

        self.set_selected_object(pp.Point2D(x_snap, y_snap), figure=self.plotter.currFig)
        if self.selected_object_temphandle is not None:
            self.selected_object_temphandle.remove()
        self.selected_object_temphandle = self.RS_lines[ev.inaxes].ax.plot(x_snap, y_snap, 'go')[0]
        self.fig.canvas.draw()
        print("Last output = (index, distance, point)")
        print("            = (%i, %.3f, (%.3f, %.3f))" % (data[0], data[1],
                                                          x_snap, y_snap))

        self.fig.canvas.mpl_disconnect(self.mouse_cid)
        self.mouse_wait_state_owner = None

    def mouse_event_user_function(self, ev):
        if ev.inaxes not in self.cb_axes:
            print('Must select axes for which callbacks have been defined.')
            return

        print("\n(%.4f, %.4f)" %(ev.xdata, ev.ydata))
        fs, fvecs = self.user_func(ev.xdata, ev.ydata)
        print(fs)
        self.last_output = (fs, fvecs)
        self.set_selected_object(pp.Point2D(ev.xdata, ev.ydata), figure=self.plotter.currFig)
        if self.selected_object_temphandle is not None:
            self.selected_object_temphandle.remove()
        self.selected_object_temphandle = ev.inaxes.plot(ev.xdata, ev.ydata, 'go')[0]
        self.fig.canvas.draw()
        self.fig.canvas.mpl_disconnect(self.mouse_cid)
        self.mouse_wait_state_owner = None

    def assign_user_func(self, func):
        self.user_func = func

    def setup(self, arrPlots, size=None, shape=None, with_times=True, basic_widgets=True):
        if shape is None:
            numrows = 1
            numcols = 1
            for k in arrPlots.keys():
                if int(k[0]) > numrows:
                    numrows = int(k[0])
                if int(k[1]) > numcols:
                    numcols = int(k[1])

            shape = [numrows, numcols]

        self.plotter.arrange_fig(shape, arrPlots)
        self.build_plotter(figsize=size, with_times=with_times, basic_widgets=basic_widgets)

    def declare_in_context(self, con_obj):
        # context_changed flag set when new objects created and unset when Generator is
        # created with the new context code included
        self.context_changed = True
        self.context_objects[con_obj.name] = con_obj

    def setup_gen(self, name_scheme):
        name = name_scheme()

        if self.context_changed:
            self.context_changed = False
            self.make_gen(self.body_pars, name)
        else:
            try:
                self.model = gui.gen_versioner.load_gen(name)
            except:
                self.make_gen(self.body_pars, name)
            else:
                self.model.set(pars=self.body_pars)

    def set_selected_object(self, selected_object, figure= None):
        """
        Set a context_object as "selected", displaying it in bold.
        """
        try:
            selected_handle = selected_object.handle
        except AttributeError:
            selected_handle = selected_object

        fig_struct, figure = self.plotter._resolve_fig(figure)
        for subplot, subplot_struct in fig_struct.arrange.items():
            for layer in subplot_struct['layers']:
                layer_struct = self.plotter._resolve_layer(figure, layer)

                ##ISSUE: box_GUI should not have markersize field. Never gets used.
                for hname, handle in layer_struct.handles.items():
                    if selected_handle is handle:
                        layer_struct.data[hname]['selected'] = True
                        layer_struct.data[hname]['linewidth'] = 2.5
                        layer_struct.data[hname]['markersize'] = 12
                    else:
                        layer_struct.data[hname]['selected'] = False
                        layer_struct.data[hname]['linewidth'] = 1
                        layer_struct.data[hname]['markersize'] = 6

                    #handle.set_linewidth(layer_struct.data[hname]['linewidth'])
                    #if isinstance(handle, mpl.lines.Line2D):
                        #handle.set_markersize(layer_struct.data[hname]['markersize'])

        self.selected_object = selected_object
        self.plotter.show(ignore_wait = True)

    def user_update_func(self):
        """
        Function overridden in user's app, called whenever navigation keys are used to move a graphics object.
        """

    def user_pick_func(self, ev):
        """
        Function overridden in user's app, called whenever an artist is picked.
        """

    def make_gen(self, pardict, name):
        """
        Empty method, which must be overridden by the user. Requires a dictionary object containing model parameters and their values (pardict)
        as well as a string (name). Uses elements of pardict to create an instance of PyDSTool.common.args,
        DSargs, which must be used to assign a model to the diagnosticGUI object at the end of the function:

        self.model = self.gen_versioner.make(DSargs)

        In order to create events in the model, a targetlang must also be assigned for PyDSTool:

        targetlang = self.gen_versioner._targetlangs[self.gen_versioner.gen_type]

        """
        raise NotImplementedError

class context_object(object):
    # Abstract base class
    pass

class data_GUI(object):
    """
    Currently, this class just makes it convenient to retrieve important fovea properties (such as the layer and name)
    of data in dstruct, given a matplotlib object (returned by on_pick).

    Right now data_GUIs are created on the spot when the selected_object changes and never referred to again.
    In the future, data_GUIs should be persistent, created when data is added to the plotter, and finally drawn in .build_layer
    using attributes stored here (not the other way around), much like how shape_GUIs are created.

    May want to make this a subclass of context_object as well.
    """
    def __init__(self, name, handle, layer, gui):
        self.name = name
        self.handle = handle
        self.layer = layer
        self.gui = gui

        self.__str__ = self.__repr__

    def getData(self):
        """
        Recovers data associated with this data_GUI from the fig_struct.
        """
        fig_struct, fig = self.gui.plotter._resolve_fig(None)

        data = fig_struct.layers[self.layer].data[self.name]['data']

        return data

    def __repr__(self):
        return "data_GUI(%s: %s in %s)" %(self.name, self.handle.__str__(), self.layer)

class domain_GUI(context_object):
    """
    Irregular polygon context objects created by domain2D. Not fully implemented.
    (e.g., missing other methods used in shape_GUI)

    ISSUE: Creating domainGUI's causes other context objects to behave strangely. For instance,
    if a line_GUI already exists, then a domain_GUI is created, the line_GUI will become bold
    and will no longer update.
    """
    def __init__(self, gui, coords, layer='gx_objects', name= None, subplot= None):
        self.name = name
        self.layer = layer
        self.handle = []
        self.gui = gui

        self.gui.declare_in_context(self)

        self.gui.plotter.add_obj(coords, mpl.lines.Line2D, name= self.name, layer= layer, style= 'y-')

        self.gui.plotter.show(ignore_wait = True)

class shape_GUI(context_object):
    def __init__(self, gui, pt1, pt2, layer='gx_objects', name= None, subplot=None):
        self.gui = gui
        self.handle = [] #mpl line or rectangle object will be stored here when draw to the axes.

        xnames = pt1.coordnames
        if pt2.coordnames != xnames:
            raise ValueError("Coordinate name mismatch")
        x1, y1 = pt1.coordarray
        x2, y2 = pt2.coordarray

        fig_struct, figure = self.gui.plotter._resolve_fig(None)
        try:
            self.gui.plotter._resolve_layer(figure, layer)
        except KeyError:
            if subplot is None:
                raise ValueError("Must specify a subplot if layer %s doesn't already exist."%layer)
            self.gui.plotter.add_layer(layer, subplot = subplot, kind = 'obj') #Set to active layer? True.
            print("Created layer %s to support Line_GUI object"%layer)

        [x1, x2, y1, y2] = self.order_points(x1, x2, y1, y2)

        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.dy = self.y2-self.y1
        self.dx = self.x2-self.x1

        self.extra_fnspecs = {}
        self.extra_pars = {}
        self.extra_auxvars = {}
        self.extra_events = []

        self.layer = layer
        if name is None:
            name = 'untitled1'

            #Increment the number if "untitled" already in use.
            while name in list(self.gui.context_objects.keys()):
                name = name[:-1]+str(int(name[-1])+1)
        else:
            if name in list(self.gui.context_objects.keys()):
                raise NameError('Name already in use by another context object.')

        self.name = name

        # declare self to GUI
        self.gui.declare_in_context(self)

    def show(self, draw= True):
        fig_struct, figure = self.gui.plotter._resolve_fig(None)
        dstruct = fig_struct.layers[self.layer]['data'][self.name]
        self.gui.plotter.set_data(self.layer, data={self.name: {'data': dstruct['data'], 'obj':dstruct['obj'],
                                                               'style':dstruct['style'], 'subplot':dstruct['subplot'],
                                                               'selected':dstruct['selected'],'linewidth':dstruct['linewidth'],
                                                               'markersize':dstruct['markersize'],'display': True}})
        if draw:
            self.gui.plotter.show(ignore_wait = True)

    def unshow(self, draw= True):
        fig_struct, figure = self.gui.plotter._resolve_fig(None)
        dstruct = fig_struct.layers[self.layer]['data'][self.name]
        self.gui.plotter.set_data(self.layer,
                                 data={self.name: {'data': dstruct['data'], 'obj':dstruct['obj'],
                                                   'style':dstruct['style'], 'subplot':dstruct['subplot'],
                                                   'selected': dstruct['selected'], 'linewidth':dstruct['linewidth'],
                                                   'markersize':dstruct['markersize'], 'display': False}})

        if draw:
            self.gui.plotter.show()

    def remove(self, draw= True):
        self.unshow(draw= draw)
        fig_struct, figure = self.gui.plotter._resolve_fig(None)
        fig_struct.layers[self.layer]['handles'].pop(self.name)
        self.gui.plotter.set_layer(self.layer, handles = fig_struct.layers[self.layer]['handles'])
        self.gui.context_objects.pop(self.name)

    def update(self, name= None, x1= None, y1= None, x2= None, y2= None):
        fig_struct, figure = self.gui.plotter._resolve_fig(None)
        show = False

        if name is not None:
            if name in list(self.gui.context_objects.keys()):
                print("Error: There already exists a context object with that name.")
                return

            for field in ['handles', 'data', 'trajs']:
                fig_struct.layers[self.layer][field][name] = fig_struct.layers[self.layer][field][self.name]
                del(fig_struct.layers[self.layer][field][self.name])

            self.gui.context_objects[name] = self.gui.context_objects.pop(self.name)
            self.name = name

        ##ISSUE: ['data'][1][0] and ['data'][1][1] represent width/height in box_GUI, NOT coords.
        ## ['data'][0][0] and ['data'][0][1] represent an x and y coord in box_GUI, NOT the x coords.
        ##Function may need clean-up due to difference in internal representation of these two objects.
        if isinstance(self, line_GUI):
            if y1 is not None:
                self.y1 = y1
                fig_struct.layers[self.layer]['data'][self.name]['data'][1][0] = y1
                self.b = self.y2 - self.x2*self.m
                show = True

            if x1 is not None:
                self.x1 = x1
                fig_struct.layers[self.layer]['data'][self.name]['data'][0][0] = x1
                self.b = self.y2 - self.x2*self.m
                show = True

            if x2 is not None:
                self.x2 = x2
                fig_struct.layers[self.layer]['data'][self.name]['data'][0][1] = x2
                self.b = self.y2 - self.x2*self.m
                show = True

            if y2 is not None:
                self.y2 = y2
                fig_struct.layers[self.layer]['data'][self.name]['data'][1][1] = y2
                self.b = self.y2 - self.x2*self.m
                show = True

            self.m = (self.y2 - self.y1)/(self.x2 - self.x1)

        if isinstance(self, box_GUI):
            if y1 is not None:
                self.y1 = y1
                fig_struct.layers[self.layer]['data'][self.name]['data'][0][1] = y1
                show = True

            if x1 is not None:
                self.x1 = x1
                fig_struct.layers[self.layer]['data'][self.name]['data'][0][0] = x1
                show = True

            if x2 is not None:
                self.dx = x2
                fig_struct.layers[self.layer]['data'][self.name]['data'][1][0] = x2
                show = True

            if y2 is not None:
                self.dy = y2
                fig_struct.layers[self.layer]['data'][self.name]['data'][1][1] = y2
                show = True

        self.gui.user_update_func()

        if show:
            self.gui.plotter.show(ignore_wait = True)

    def make_event_def(self, uniquename, dircode=0):
        """
        make_event_def was created for line_GUIs, before box_GUI was a distinct object. This method has not
        been tested for box_GUIs and will probably behave unexpectedly.

        make_event_def should probably also be called in .update, if an event already exists. Right now,
        if a line is created, then moved, the event will behave as though there were a line in the original spot.
        """
        fig_struct, figure = self.gui.plotter._resolve_fig(None)

        #for field in ['handles', 'data', 'trajs']:
            #fig_struct.layers[self.layer][field][uniquename] = \
                #fig_struct.layers[self.layer][field].pop(self.name)

        ##Update changes the diagnosticGUI object.
        #self.gui.plotter.figs[figure] = fig_struct

        #self.name = uniquename
        res = pp.make_distance_to_line_auxfn('exit_line_'+uniquename, 'exit_fn_'+uniquename,
                                          ['x', 'y'], True)

        parname_base = 'exit_line_%s_' %uniquename
        self.extra_pars[parname_base+'p_x'] = self.x1
        self.extra_pars[parname_base+'p_y'] = self.y1
        self.extra_pars[parname_base+'dp_x'] = self.dx
        self.extra_pars[parname_base+'dp_y'] = self.dy
        self.extra_fnspecs.update(res['auxfn'])
        targetlang = \
            self.gui.gen_versioner._targetlangs[self.gui.gen_versioner.gen_type]
        self.extra_events = [dst.Events.makeZeroCrossEvent(
                                              expr='exit_fn_%s(x,y)' %uniquename,
                                              dircode=dircode,
                                              argDict={'name': 'exit_ev_%s' %uniquename,
                                                       'eventtol': 1e-8,
                                                       'eventdelay': 1e-3,
                                                       'starttime': 0,
                                                       'precise': True,
                                                       'active': True,
                                                       'term': False},
                                              varnames=('x', 'y'),
                                              fnspecs=res['auxfn'],
                                              parnames=res['pars'],
                                              targetlang=targetlang
                                              )]


class box_GUI(shape_GUI):
    """
    Box of interest context_object for GUI
    """
    def __init__(self, gui, pt1, pt2, layer='gx_objects', subplot=None, name=None, select= True):

        shape_GUI.__init__(self, gui, pt1, pt2, layer='gx_objects', subplot=subplot, name= name)

        #if self.x1 > self.y2:
            #xt = self.x1
            #self.x1 = self.x2
            #self.x2 = xt

        #if self.y1 > self.y2:
            #yt = self.y1
            #self.y1 = self.y2
            #self.y2 = yt

        if select:
            self.gui.set_selected_object(self, figure= self.gui.plotter.currFig)
            print("Created box and moved to currently selected object")
        self.gui.plotter.add_obj(np.array([[self.x1, self.y1],[self.dx, self.dy]]), mpl.patches.Rectangle,
                                layer=layer, subplot=subplot, style= None, name=self.name, force= True, display= True)

        self.show()

    def __repr__(self):
        if self.extra_events == []:
            ev_str = '(no event)'
        else:
            ev_str = '(with event)'
        return "box_GUI(%.3f, %.3f, %.3f, %.3f) - '%s' %s" %(self.x1, self.y1, \
                                          self.x2, self.y2, self.name, ev_str)

    def pin_contents(self, traj, coorddict):
        """
        Determine if a trajectory passes through the box object and add that segment of the trajectory
        to a layer. Uses coorddict format of add_data_points.

        Current version only takes horizontal slices of the trajectory, preserving all change in the y-direction.
        Should be adapted to 2D. It would be preferable if traj was not required as a param, and found with
        a test for containment.

        Note: This method may be redundant when better event creation/handling has been implemented
        for box_GUIs.
        """
        for var, params in coorddict.items():
            pts = traj.sample(tlo= self.x1, thi= self.x1 + self.dx)
            pts[params['x']] = pts[params['x']] - self.x1

            xs = pts[params['x']]
            ys = pts[var]

            self.gui.plotter.add_data([xs, ys], layer= params['layer'],
                                 name= params['name'], style = params['style'], traj= pts, force= True)

    def order_points(self, x1, x2, y1, y2):
        if x1 > x2:
            xt = x1
            x1 = x2
            x2 = xt

        if y1 > y2:
            yt = y1
            y1 = y2
            y2 = yt

        return [x1, x2, y1, y2]


class line_GUI(shape_GUI):
    """
    Line of interest context_object for GUI
    """
    def __init__(self, gui, pt1, pt2, layer='gx_objects', subplot=None, name= None, select= True):

        shape_GUI.__init__(self, gui, pt1, pt2, layer='gx_objects', subplot=subplot, name= name)

        self.length = np.linalg.norm((self.dx, self.dy))
        # angle relative to horizontal, in radians
        self.ang = atan2(self.dy,self.dx)
        self.ang_deg = 180*self.ang/pi

        # slope and y intercept
        self.m = (self.y2 - self.y1)/(self.x2 - self.x1)
        self.b = self.y2 - self.x2*self.m

        if select:
            self.gui.set_selected_object(self, figure= self.gui.plotter.currFig)
            print("Created line and moved to currently selected object")

        # actual MPL line object handle
        self.gui.plotter.add_obj(np.array([[self.x1, self.x2],[self.y1, self.y2]]), mpl.lines.Line2D, layer=layer, subplot=subplot,
                           style= None, name=self.name, force= True, display= True)
        self.show()

    def __repr__(self):
        if self.extra_events == []:
            ev_str = '(no event)'
        else:
            ev_str = '(with event)'
        return "line_GUI(%.3f, %.3f, %.3f, %.3f) - '%s' %s" %(self.x1, self.y1, \
                                          self.x2, self.y2, self.name, ev_str)

    def distance_to_pos(self, dist):
        """
        Calculate absolute (x,y) position of distance dist from (x1,y1) along line
        """
        return self.fraction_to_pos(self, dist/self.length)


    def fraction_to_pos(self, fraction):
        """
        Calculate absolute (x,y) position of fractional distance (0-1) from (x1,y1) along line
        """
        return np.array([self.x1+fraction*self.dx, self.y1+fraction*self.dy])

    def order_points(self, x1, x2, y1, y2):
        if x1 > x2:
            # ensure correct ordering for angles
            xt = x1
            x1 = x2
            x2 = xt
            yt = y1
            y1 = y2
            y2 = yt

        return [x1, x2, y1, y2]

    def points(self):
        """
        Return 2D list of points along the line.
        """
        m = (self.y2 - self.y1)/(self.x2 - self.x1)
        b = self.y2 - self.x2*m
        xs = list(range(int(self.x1), int(self.x2))) #Casting as ints may be risky.
        #ys = map(lambda x: x*m + b, xs)
        ys = [x*m+b for x in xs]

        return [xs, ys]


class tracker_GUI(object):
    """
    Abstract base class
    """
    pass


class tracker_textconsole(tracker_GUI):
    """
    Auto-updating text consoles that are connected to a diagnostic GUI.
    """
    def __init__(self):
        self.figs = {}
        self.sim = None
        self.calc_context = None

    def __call__(self, calc_context, fignum, attribute_name):
        self.sim = calc_context.sim
        self.calc_context = calc_context
        old_toolbar = plt.rcParams['toolbar']
        plt.rcParams['toolbar'] = 'None'
        fig = plt.figure(fignum, figsize=(2,6)) #, frameon=False)
        plt.rcParams['toolbar'] = old_toolbar
        if fignum in self.figs:
            self.figs[fignum].tracked.append(attribute_name)
        else:
            self.figs[fignum] = args(figure=fig, tracked=[attribute_name])
        self.sim.tracked_objects.append(self)

    def show(self):
        for fignum, figdata in self.figs.items():
            fig = plt.figure(fignum)
            ax = plt.axes([0., 0., 1., 1.], frameon=False, xticks=[],yticks=[])
            #figdata.figure.clf()
            ax.cla()
            ax.set_frame_on(False)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            wspace = self.calc_context.workspace
            for tracked_attr in figdata.tracked:
                for i, (obj_name, obj) in enumerate(wspace.__dict__.items()):
                    if obj_name[0] == '_':
                        # internal name, ignore
                        continue
                    try:
                        data = getattr(obj, tracked_attr)
                    except Exception as e:
                        print("No attribute: '%s' in object in workspace '%s'" % (tracked_attr, wspace._name))
                        raise
                    plt.text(0.05, 0.05+i*0.04, '%s: %s = %.4g' % (obj_name, tracked_attr, data))
            plt.title('%s measures of %s (workspace: %s)'%(self.calc_context.sim.name, tracked_attr,
                                                           _escape_underscore(self.calc_context.workspace._name)))
            fig.canvas.set_window_title("Fig %i, Workspace %s" % (fignum, self.calc_context.workspace._name))
        #plt.show()



class tracker_plotter(tracker_GUI):
    """
    Auto-updating plots that are connected to a diagnostic GUI.
    """
    def __init__(self):
        self.figs = {}
        self.sim = None
        self.calc_context = None

    def __call__(self, calc_context, fignum, xstr, ystr, style):
        self.sim = calc_context.sim
        self.calc_context = calc_context
        fig = plt.figure(fignum)
        new_track = args(xstr=xstr, ystr=ystr, style=style)
        if fignum in self.figs:
            self.figs[fignum].tracked.append(new_track)
        else:
            self.figs[fignum] = args(figure=fig, tracked=[new_track])
        self.sim.tracked_objects.append(self)

    def show(self):
        for fignum, figdata in self.figs.items():
            fig = plt.figure(fignum)
            ax = plt.gca()
            #figdata.figure.clf()
            ax.cla()
            wspace = self.calc_context.workspace
            for tracked in figdata.tracked:
                try:
                    xdata = getattr(wspace, tracked.xstr)
                except Exception as e:
                    print("Failed to evaluate: '%s' in workspace '%s'" % (tracked.xstr, wspace._name))
                    raise
                try:
                    ydata = getattr(self.calc_context.workspace, tracked.ystr)
                except Exception as e:
                    print("Failed to evaluate: '%s' in workspace '%s'" % (tracked.ystr, wspace._name))
                    raise
                plt.plot(xdata, ydata,
                     tracked.style, label=_escape_underscore(tracked.ystr))
            plt.legend()
            plt.title('%s measures vs %s (workspace: %s)'%(self.calc_context.sim.name, tracked.xstr,
                                                           _escape_underscore(self.calc_context.workspace._name)))
            fig.canvas.set_window_title("Fig %i, Workspace %s" % (fignum, self.calc_context.workspace._name))
        #plt.show()

class tracker_manager(object):
    """
    Track different quantities from different calc contexts in different
    figures. Cannot re-use same figure with different contexts.
    """
    def __init__(self):
        self.tracked = {}
        # currently, all_figs does not automatically release figures if
        # contexts are deleted or replaced
        self.all_figs = []

    def __call__(self, calc_context, fignum, plot_metadata=None,
                 text_metadata=None):
        """
        plot_metadata (default None) = (xstr, ystr, style)
        *or*
        text_metadata (default None) = attribute_name

        If text_metadata used, the tracker object is assumed to be a
        textconsole type that accesses declared python objects to access an
        attribute

        """
        valid_input = plot_metadata is None or text_metadata is None
        if not valid_input:
            raise ValueError("Only use one of plot or text metadata arguments")
        try:
            xstr, ystr, style = plot_metadata
        except TypeError:
            # None is not iterable
            track_plot = False
            attribute_name = text_metadata
        else:
            track_plot = True
        if calc_context in self.tracked:
            if track_plot:
                # tracker_plotter type
                self.tracked[calc_context](calc_context, fignum, xstr, ystr, style)
            else:
                # tracker_textconsole type
                self.tracked[calc_context](calc_context, fignum, attribute_name)
            if fignum not in self.all_figs:
                self.all_figs.append(fignum)
        else:
            if fignum in self.all_figs:
                raise ValueError("Figure number %i already in use" % fignum)
            if track_plot:
                tp = tracker_plotter()
                tp(calc_context, fignum, xstr, ystr, style)
            else:
                tp = tracker_textconsole()
                tp(calc_context, fignum, attribute_name)
            self.tracked[calc_context] = tp


    def show(self):
        for tp in self.tracked.values():
            tp.show()


##def track_attribute(attr_name):
##    """
##    Create decorator to track named attribute used in a calculation
##    """
##    def decorator(fn):
##        obj =
##        tracker.track_list = [getattr(obj, attr_name)]
##        return fn
##    return decorator

def _escape_underscore(text):
    """
    Internal utility to escape any TeX-related underscore ('_') characters in mpl strings
    """
    return text.replace('_', '\_')

# ---------------------------------------------------------

#global gui, tracker

# singleton pattern

# plotter is not a globally visible object outside of this module
plotter = Plotter()

gui = diagnosticGUI(plotter)

tracker = tracker_manager()
