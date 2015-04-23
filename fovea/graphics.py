"""
Graphical user interface and plotting tools for dynamical systems

Rob Clewley, 2015
based on work by:
Bryce Chung and Rob Clewley, 2012


==== Usage notes:
Plotting styles can be given either as a string or as a dictionary of
  style kwargs suitable for plot command

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RectangleSelector
import numpy as np
from copy import copy
from math import *
import hashlib, time
import euclid as euc

from PyDSTool import args, numeric_to_traj, Point
import PyDSTool.Toolbox.phaseplane as pp
from PyDSTool.Toolbox.phaseplane import Point2D
# for potentially generalizable functions and classes to use
import PyDSTool as dst

# local imports
from common import *
import domain2D as dom


# ----------------------------------------------


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
        return Point2D(p1), Point2D(p2)
    else:
        raise ValueError("No intersection")



class plotter2D(object):

    colors = ['b', 'g', 'r', 'c', 'm', 'k', 'y']

    def __init__(self):
        self.clean()

    def clean(self):
        """
        Delete all figure data (doesn't clear figures)
        """
        self.figs = {}
        self._max_fig_num = 0
        self.active_layer = None
        self.currFig = None
        # record whether ever called show()
        self.shown = False

    def auto_scale_domain(self, figure=None):
        """
        Set domain limits to that of the data in all layers
        with the greatest extent.
        """
        # ISSUE: figure not used! It should select which
        # data is used below, and the figure object should
        # be selected to impose the extent on

        # initial values
        x_extent = [0,0]
        y_extent = [0,0]

        found_fig = False
        for figName, fig in self.figs.items():
            if figure != figName:
                continue
            else:
                found_fig = True
            for layerName, layer in fig.layers.items():
                for dName, d in layer['data'].items():
                    x_extent[0] = min(min(d[0][0]), x_extent[0])
                    x_extent[1] = max(max(d[0][0]), x_extent[1])
                    y_extent[0] = min(min(d[0][1]), y_extent[0])
                    y_extent[1] = max(max(d[0][1]), y_extent[1])

        if not found_fig:
            raise ValueError("No such figure")
        fig.domain = (x_extent, y_extent)
        plt.figure(fig.fignum)
        plt.xlim(x_extent)
        plt.ylim(y_extent)

    def set_active_layer(self, layer, figure=None):
        """
        Sets the active_layer attribute to be the named layer struct
        in the (optionally) given figure (defaults to Master)
        """
        fig_struct, figure = plotter._resolveFig(figure)
        self.active_layer = fig_struct.layers[layer]


    # ISSUE: This method is never called!
##    def build(self, wait=False, figure=None, labels=None, autoscale=True):
##        """
##        Separate function for building figures to optimize simulation speeds.
##        Call build() when user wants figure(s) to be displayed.
##
##        wait        Will automatically pause for user command after plotting figures when set to True; default False
##        figure      Allows user to specify which figure to plot; overwrites figure.display settings
##        labels      Specify whether to display labels on the plots and at which level
##                    None | layers | data
##        """
##
##        for figName in self.figs:
##
##            fig_struct = self.figs[figName]
##
##            if not fig_struct.display:
##                continue
##
##            fig = plt.figure(fig_struct.fignum)
##            #plt.title(fig_struct.title)
##
##            if len(fig_struct.arrange) == 0:
##                ax = fig.add_subplot(1,1,1)
##
##                for layName in fig_struct.layers:
##                    self._buildLayer_(figName, layName, ax)[0]
##                ax.set_xlabel(fig_struct.xlabel)
##                ax.set_ylabel(fig_struct.ylabel)
##                ax.set_xlim(fig_struct.xdom)
##                ax.set_ylim(fig_struct.ydom)
##                ax.autoscale(enable=autoscale)
##            else:
##                shape = fig_struct.shape
##                for ix, subplot in fig_struct.arrange.items():
##                    ax = fig.add_subplot(shape[0], shape[1], shape[1]*(int(ix[0])-1)+int(ix[1]))
##                    try:
##                        axesLabels = subplot['axes']
##                        if axesLabels is not None:
##                            ax.set_xlabel(axesLabels[0])
##                            ax.set_ylabel(axesLabels[1])
##                    except:
##                        raise ValueError("Error setting axis label for subplot")
##                    try:
##                        plt.title(subplot['name'])
##                    except:  #### What exception is caught here?
##                        pass
##                    layer = subplot['layers']
##                    if isinstance(layer, list):
##                        if len(layer) > 0:
##                            for layName in layer:
##                                self._buildLayer_(figName, layName, ax)
##                        else:
##                            continue
##                    else:
##                        self._buildLayer_(figName, layer, ax)
##                    plt.autoscale(enable=autoscale)
##
##            fig.canvas.draw()
##            self.figs[figName].window = fig
##
##        if wait:
##            pause = raw_input("Press <Enter> to continue...")


    def show_legends(self, figure=None, subplot=None):
        """
        Show all figure legends of visible data layers to stdout.
        Option to filter results by sub-plot string code or name, e.g. '21'
        """
        fig_struct, figure = self._resolveFig(figure)
        arrPlots = fig_struct.arrange

        for ixstr, spec in arrPlots.items():
            if subplot is not None and subplot not in \
                              [ixstr, spec['name']]:
                continue
            layer_info = spec['layers']
            for layer in layer_info:
                lay = fig_struct.layers[layer]
                if not lay.display or lay.kind != 'data':
                    continue
                print("\nLAYER: %s" % str(layer))
                print("  sub-plot: %s" % ixstr)
                print("  style %s" % lay.style)
                print("  axes: %s - %s" % (lay.axes_vars[0], lay.axes_vars[1]))
                print("  data:")
                for dname, dstruct in lay.data.items():
                    if not dstruct['display']:
                        continue
                    print("     name: %s, style: %s" % (dname, dstruct['style']))



    ## Figure Management Tools ##

    # ISSUE: make compound method names with _ not camelCase
    def addFig(self, label, title="", xlabel="", ylabel="", tdom=None,
               domain=None, display=True):
        """
        User can add figures to plotter for data shown in different windows

        label        String to specify name of figure
        title        String to specify title to display on figure
        xlabel       String to specify label on x-axis
        ylabel       String to specify label on y-axis
        tdom         time domain (if applicable, and optional)
        display      Setting to determine whether or not figure is plotted when
                       plotter is built (defaults to True)
        """
        # Check to see if label already exists
        if self.figs.has_key(label):
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


    def copyFig(self, newFig, oldFig):
        """
        Duplicates all figure, layer, and data information without arrangement
        information.
        """
        # Check to see if label already exists #
        if self.figs.has_key(newFig):
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
            self.addLayer(layer, display=oldLay.display, zindex=oldLay.zindex, style=oldLay.style)


    def setFig(self, label=None, **kwargs):
        """
        Sets current figure given by parameter 'label', if given.

        Also alows user to set properties of the named (or current) figure.
        See addFig() for properties to set
        """
        # Check to see that figure exists #
        if label == None and self.currFig != None:
            label = self.currFig
        elif label == None and self.currFig == None:
            raise ValueError("Must set current figure")
        elif label != None:
            if not self.figs.has_key(label):
                raise KeyError("Figure does not exist with specified label!")

        self.currFig = label

        for kw in kwargs:
            # Check to see that keyword is valid #
            if not self.figs[label].has_key(kw):
                raise KeyError("Unable to set figure property: Invalid keyword argument")

            self.figs[label][kw] = kwargs[kw]


    def clearFig(self, label):
        """
        Ignore if no figure with given label is found.
        """
        if self.figs.has_key(label):
            for layer in self.figs[label].layers:
                self.figs[label].layers[layer].data = {}

            if self.currFig == label:
                self.currFig = None


    ## Layer Management Tools ##

    def arrangeFig(self, shape, arrPlots, figure=None):
        """
        Separate layers in a figure into different subplots.

        shape      [rows,cols] where rows is number of rows starting from 1
                   and cols is number of cols starting from 1
        arrPlots   Dict of dicts of subplots indexed by row-col.  Each
                   subplot dict has the following keys:

                   name        Name of the subplot to appear in figure
                   layers      List of layers to be displayed in subplot
                   axes        Name of axes, [xlabel, ylabel]
                   scale       Scale of axes [(xlo, xhi), (ylo, yhi)]
                                (will NOT overwrite any already declared for a layer)

        figure     Figure name to apply changes to. If left blank, current
                   figure is used

        """
        fig_struct, figure = self._resolveFig(figure)

        if len(shape) != 2:
            raise ValueError("shape must be (rows,cols)")

        # make copy of arrPlots in case change singleton layer names to lists
        arrPlots = copy(arrPlots)

        for ixstr, spec in arrPlots.items():
            if int(ixstr[0])*int(ixstr[1]) > shape[0]*shape[1]:
                raise ValueError("Position does not exist in subplot arrangement.")

            layer_info = spec['layers']
            if len(spec['axes_vars']) > 2:
                raise ValueError("Cannot have more than two axis titles.")
            else:
                axes_vars = spec['axes_vars']

            if not isinstance(layer_info, list):
                # will be a singleton string layer name
                layer_info = [layer_info]
                # resave this into the version to be stored
                spec['layers'] = layer_info

            for layName in layer_info:
                if layName not in fig_struct.layers.keys():
                    raise KeyError("Layer not created in figure: %s" % str(layName))
                fig_struct.layers[layName].axes_vars = axes_vars
                if 'scale' in spec and fig_struct.layers[layName].scale is None:
                    # i.e. don't overwrite an existing, specific layer's scale
                    fig_struct.layers[layName].scale = checked_scale(spec['scale'])

        fig_struct.shape = shape
        fig_struct.arrange = arrPlots


    def addLayer(self, layer_name, figure=None, **kwargs):
        """
        User method to add data sets to a layer in a figure

        figure  name
        layer   name
        display toggle Boolean
        kind    a user-defined kind, e.g. 'data', 'vline', 'hline', 'epoch', etc.
        scale   a pair of axes scale pairs or None
        zindex  (currently unused)
        """
        fig_struct, figure = self._resolveFig(figure)

        # Check to see layer does not already exist
        if fig_struct.layers.has_key(layer_name):
            raise KeyError("Layer name already exists in figure!")

        # default plot style generated from list of colors
        color = self.colors[len(fig_struct.layers)%len(self.colors)]
        line = '-'
        style = color+line

        layAttrs = args()
        layAttrs.data = {}
        layAttrs.display = True
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

        for kw in kwargs:
            # Check to see that parameter exists in layers
            ## Possibly change to account for different properties of specific artist objects?
            if not layAttrs.has_key(kw):
                raise KeyError("Parameter is not a property of the layer.")

            layAttrs[kw] = kwargs[kw]

        fig_struct.layers[layer_name] = layAttrs


    def setLayer(self, label, figure=None, **kwargs):
        """
        Arrange data sets in a figure's layer
        """
        # figure will be the same before and after unless figure was
        # None, in which case defaults to name of master figure
        fig_struct, figure = self._resolveFig(figure)

        # Check to see that layer exists
        if not fig_struct.layers.has_key(label):
            raise KeyError("Layer does not exist in figure!")

        for kw in kwargs:
            # Check to see that parameter exists in layers
            ## Possibly change to account for different properties of specific artist objects?
            if not fig_struct.layers[label].has_key(kw):
                raise KeyError("Parameter is not a property of the layer.")

            fig_struct.layers[label][kw] = kwargs[kw]


    def addData(self, data, figure=None, layer=None, style=None, name=None, disp=True,
                force=False):
        """
        User tool to add data to a layer.
        Use force option only if known that existing data must be overwritten.
        """

        # Check to see that data is a list or array
        try:
            size = np.shape(data)
        except:
            raise TypeError("Data must be castable to a numpy array")

        # Check to see that there is an x- and y- dataset
        if size[0] != 2:
            raise ValueError("Data must contain 2 lists of data points")

        fig_struct, figure = self._resolveFig(figure)

        if layer is None:
            layer = "layer"+str(len(fig_struct.layers)+1)
            self.addLayer(figure, layer)
        elif not fig_struct.layers.has_key(layer):
            raise KeyError("Layer has not been created")

        # inherit default style from layer if not given
        if style is None:
            style = fig_struct.layers[layer].style

        # d is a dictionary mapping 'name' to a dictionary of fields for the
        # numerical plot data (key 'data'), style (key 'style'), and display
        # boolean (key 'display').
        #
        # The numerical data for d is given by the method argument also called data
        d = fig_struct.layers[layer].data

        # Check to see if data name already exists
        if d.has_key(name) and not force:
            raise KeyError("Data name already exists in layer: %s" %name)

        # Create name if none is provided
        if name is None:
            i = 1
            flag = True
            while flag:
                strName = str(layer)+str(i)
                if d.has_key(strName):
                    i = i+1
                else:
                    flag = False

            name = strName

        #print("Adding: %s >> %s" % (figure, name))
        d.update({name: {'data': data, 'style': style, 'display': disp}})
        self._updateTraj(layer, figure)


    def setPoint(self, name, pt, layer, figure=None):
        """
        Set point coordinates in given layer, specified by pt.

        Figure defaults to currFig if not specified.
        """
        fig_struct, figure = self._resolveFig(figure)

        pt_struct = fig_struct.layers[layer].data[name]
        pt_struct['data'] = [[pt.x], [pt.y]]


    def setData(self, layer, figure=None, **kwargs):
        """
        Set data properties in given layer, specified by
        the keys of the keyword arguments.

        Figure defaults to currFig if not specified.
        """
        fig_struct, figure = self._resolveFig(figure)

        # Check to see that layer exists
        if not fig_struct.layers.has_key(layer):
            raise KeyError("Layer does not exist in figure!")

        for kw in kwargs:
            # Possibly change to account for different properties of specific artist objects?
            try:
                # for dict updates, e.g. data={<line_name>: {'data': [<array_data>]}}
                fig_struct.layers[layer][kw].update(kwargs[kw]) # = kwargs[kw]
            except AttributeError:
                # for non-dict updates, e.g. display=True
                fig_struct.layers[layer][kw] = kwargs[kw]
            except KeyError:
                raise KeyError("Parameter '%s' is not a property of the layer." %kw)

        self._updateTraj(layer, figure)



    def _updateTraj(self, layer_name, figure_name):
        fig_struct = self.figs[figure_name]
        layer = fig_struct.layers[layer_name]
        #dstruct = layer.data  ???
        for name, dstruct in layer.data.items():
            # catch for when the data is individual points
            try:
                layer.trajs[name] = numeric_to_traj([dstruct['data'][1]],
                                        layer_name+'.'+name,
                                        coordnames=['y'],
                                        indepvar=dstruct['data'][0],
                                        discrete=False)
            except ValueError:
                pass


    def appendData(self, data, layer, name, figure=None):
        """
        Doc string?
        """
        fig_struct, figure = self._resolveFig(figure)

        try:
            lay = fig_struct.layers[layer]
        except KeyError:
            raise ValueError("Layer %s does not exist figure"%layer)

        try:
            dataset = lay.data[name]
        except KeyError:
            raise ValueError("Dataset %s does not exist in layer" % name)

        if isinstance(data, PyDSTool.Toolbox.phaseplane.Point2D):
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

        self._updateTraj(layer, figure)


    def addLineByPoints(self, pts, figure=None, layer=None, style=None, name=None):
        """
        Add line based on two given Point2D points
        """
        pt1, pt2 = pts
        try:
            self.addData([[pt1.x, pt2.x], [pt1.y, pt2.y]], figure=figure, layer=layer, style=style, name=name)
        except AttributeError:
            raise ValueError("First argument must be [Point2D, Point2D]")


    def addPoint(self, pt, figure=None, layer=None, style=None, name=None):
        """
        Add single Point2D point
        """
        self.addData([[pt.x], [pt.y]], figure=figure, layer=layer, style=style,
                     name=name)


    def addVLine(self, x, figure=None, layer=None, style=None, name=None):
        """
        Add vertical line
        """
        fig_struct, figure = self._resolveFig(figure)

        if fig_struct.layers[layer].scale[1] is None:
            # RC: Not good to assume [0, 1]
            self.addData([[x, x], [0, 1]], figure=figure, layer=layer, style=style, name=name)
            fig_struct.layers[layer].kind = 'vline'
        else:
            self.addData([[x, x], fig_struct.layers[layer].scale[1]], figure=figure, layer=layer, style=style, name=name)
            fig_struct.layers[layer].kind = 'vline'


    def addHLine(self, y, figure=None, layer=None, style=None, name=None):
        """
        Add horizontal line.

        NB. This should somehow be changed to use ax.axhline, which automatically
        always spans the x axis with the default coords settings.
        """
        fig_struct, figure = self._resolveFig(figure)

        sc = fig_struct.layers[layer].scale
        if sc is None:
            if self.figs[figure].tdom is None:
                tdom = sc[0]
            else:
                tdom = self.figs[figure].tdom
        else:
            if sc[0] is None:
                # try defaulting to figure tdom
                if self.figs[figure].tdom is None:
                    # RC: Not good to assume [0, 100]
                    tdom = [0, 100]
                else:
                    tdom = self.figs[figure].tdom
            else:
                tdom = sc[0]
        self.addData([tdom, [y, y]], figure=figure, layer=layer, style=style, name=name)
        fig_struct.layers[layer].kind = 'hline'

    def show(self):
        """
        Apply all figures' domain limits and refresh
        """
        for figName, fig in self.figs.items():
            f = plt.figure(fig.fignum)
            ax = f.gca()
            xdom, ydom = fig.domain
            ax.set_xlim(xdom)
            ax.set_ylim(ydom)
        plt.draw()
        if self.shown:
            plt.show()
            self.shown = True


    def updateDynamic(self, time, dynamicFns, hard_reset=False):
        """Dynamic callback functions always accept time as first argument.
        Optional second argument is hard_reset Boolean.
        """
        for figName in self.figs:
            for layer in self.figs[figName].layers:
                argsLayer = self.figs[figName].layers[layer]

                if not argsLayer.dynamic:
                    continue
                else:
                    assert isinstance(dynamicFns, dict), "Dynamic functions must be a dictionary of layer-function pairs"

                    if dynamicFns.has_key(layer):
                        #print("updateDynamic calling function: %s" % str(dynamicFns[layer]))
                        dynamicFns[layer](time, hard_reset)


    def _resolveFig(self, figure):
        """Internal utility to return a figure structure and figure name,
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


    def _buildLayer_(self, figure_name, layer_name, ax, rescale=None):
        """
        Consolidates layer information into matplotlib.artist objects
        rescale (pair of pairs) may be set if the axes' current scalings
           should be preserved, overriding the layer's set scale, if any
        """
        fig_struct = self.figs[figure_name]
        lay = fig_struct.layers[layer_name]
        lines = []

        for dkey, dstruct in lay.data.items():
            # we should have a way to know if the layer contains points
            # that may or may not already be updated *in place* and therefore
            # don't need to be redrawn
            if not dstruct['display']:
                continue

            # process user-defined style
            s = dstruct['style']

            if isinstance(s, str):
                style_as_string = True
            elif isinstance(s, dict):
                style_as_string = False

            if s == "":
                # default to black lines
                s = 'k-'

            # in case in future we wish to have option to reverse axes
            ix0, ix1 = 0, 1

            if style_as_string:
                lines.append(ax.plot(dstruct['data'][ix0], dstruct['data'][ix1], s))
            else:
                lines.append(ax.plot(dstruct['data'][ix0], dstruct['data'][ix1], **s))

            lay.lines = lines

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

        assert isinstance(objPlotter, plotter2D), \
               "plotter2D_GUI must be instantiated with a plotter2D object."

        # easy user access to these basic attributes: current time and
        # index into self.points
        self.t = None
        self.ix = None

        self.points = None
        # add Trajectory object too?

        # times is an array from trajectory points
        self.times = None

        # clipboard for a point selected in a dynamic sub-plot (e.g., phaseplane)
        self.clipboardPt = None

        # Captured points in all sub-plots at current time, from GUI button
        self.capturedPts = {}

        # internal verbosity level
        self.verbose = verbose_level

        if points is not None:
            self.addPoints(points)

        # ---------------
        # Internal stuff
        # ---------------

        # masterWin is the figure handle for main GUI window
        self.masterWin = None

        # lookup dictionary mapping axes handles to layer name lists
        # and sub-plot index strings
        self.subplot_lookup = {}

        self.timePlots = []
        self.timeLines = []

        self.plotter = objPlotter
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


    def addDataPoints(self, points):
        """User provides the pointset for the data to be investigated
        """
        self.points = points
        self.times = points['t']


    def buildPlotter2D(self, figsize=None, with_times=True):
        """
        Create time bar widget.
        Create capture points widget.
        Also closes all figure windows and resets masterWin attribute.

        Optional size is a pair of figure screen size measurements in inches
        """

        plt.close('all')
        self.masterWin = None

        for figName, fig_struct in self.plotter.figs.items():
            if figsize is not None:
                fig_handle = plt.figure(fig_struct.fignum, figsize=figsize)
                # fig sizing later, with fig.set_figsize_inches doesn't seem
                # to work properly unless put in the initialization call
            else:
                fig_handle = plt.figure(fig_struct.fignum)
            fig_handle.canvas.set_window_title(fig_struct.title + " : Master window")
            if figName != 'Master':
                continue
            ##### Set up master window controls
            plt.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.1,
                               wspace=0.2, hspace=0.23)

            self.masterWin = fig_handle

            ## Time bar controls time lines in figures
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

            ## Capture point button in lower left
            captureButton = Button(plt.axes([0.055, 0.02, 0.08, 0.03]), 'Capture Point')
            self.widgets['capturePoint'] = captureButton

            ## Refresh button
            refreshButton = Button(plt.axes([0.005, 0.02, 0.045, 0.03]), 'Refresh')
            self.widgets['refresh'] = refreshButton

            ## Go back to last point button
            backButton = Button(plt.axes([0.005, 0.06, 0.045, 0.03]), 'Back')
            self.widgets['goBack'] = backButton

            ## Build up each subplot, left to right, top to bottom
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
                        # singleton string layer name
                        layer_info = [layer_info]

                    # ISSUE: 'scale' should probably be 'domain' or 'extent'
                    # and the titling and labeling should happen in plotter2D
                    try:
                        scale = subplot_struct['scale']
                    except KeyError:
                        subplot_struct['scale'] = None
                        scale = None
                    ax = fig_handle.add_subplot(shape[0], shape[1], shape[1]*i + j+1)
                    self.subplot_lookup[ax] = (figName, layer_info, ixstr)
                    subplot_struct['axes_obj'] = ax

                    ax.set_title(subplot_struct['name'])
                    axes_vars = subplot_struct['axes_vars']
                    ax.set_xlabel(axes_vars[0])
                    ax.set_ylabel(axes_vars[1])

                    if len(layer_info) > 0:
                        for layName in layer_info:
                            self.plotter._buildLayer_(figName, layName, ax)
                            if self.dynamicPlotFns.has_key(layName):
                                self.dynamicPlots[layName] = ax
                                # initialize the layer's dynamic stuff
                                self.dynamicPlotFns[layName](self.t)
##                    else:
##                        ax.plot([0], [0])   # why?

                    # add vertical time line in all time plots
                    # (provide option for user to specify as time or t)
                    if axes_vars[0].lower() in ["time", 't']:
                        self.timeLines.append(ax.axvline(x=self.t, color='r',
                                                         linewidth=3, linestyle='--'))
                        self.timePlots.extend(layer_info)

                    # ISSUE: these should be built into plotter2D's figure domains instead
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

            # ISSUE: This should probably be moved to plotter2D
            fig_handle.canvas.draw()

        # activate callbacks
        self.widgets['capturePoint'].on_clicked(self.capturePoint)
        self.widgets['refresh'].on_clicked(self.refresh)
        self.widgets['goBack'].on_clicked(self.goBack)
        if with_times:
            self.widgets['timeBar'].on_changed(self.updatePlots)
            self.widgets['minus_dt'].on_clicked(self.minus_dt)
            self.widgets['plus_dt'].on_clicked(self.plus_dt)

        # activate general mouse click callbacks
        evMouseDown = fig_handle.canvas.mpl_connect('button_press_event', self.mouseDownFn)
        evMouseUp = fig_handle.canvas.mpl_connect('button_release_event', self.mouseUpFn)
        evMouseMove = fig_handle.canvas.mpl_connect('motion_notify_event', self.mouseMoveFn)
        evKeyOn = fig_handle.canvas.mpl_connect('key_press_event', self.modifier_key_on)
        evKeyOff = fig_handle.canvas.mpl_connect('key_release_event', self.modifier_key_off)


    def buildLayers(self, layer_list, ax, rescale, figure=None):
        """
        Convenience function to group layer refresh/build calls
        """
        fig_struct, figure = self.plotter._resolveFig(figure)
        for layer_name in layer_list:
            self.plotter._buildLayer_(figure, layer_name, ax, rescale)


    def clearAxes(self, subplot, figure=None):
        """
        Clears lines and points in sub-plot axes (given as an axis
        object or sub-plot string name) without removing
        other items such as title, labels.
        """
        fig_struct, figure = self.plotter._resolveFig(figure)
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
                ixstr = self.subplot_lookup[subplot][2]
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


    # standard callback functions for the GUI display

    def mouseDownFn(self, ev):
        self._mouseUp = False
        #print "mouseDownFn", self._mouseUp

    def mouseUpFn(self, ev):
        # NOTE: zoom dragging will not get completed before this callback
        # so trying to refresh as a result of zoom here will fail
        self._mouseUp = True
        #print "mouseUpFn", self._mouseUp
        # if in a time-based sub-plot and not dragging
        do_get = False
        if not self._mouseDrag:
            do_get = True
            # check widget axes
            for w in self.widgets.values():
                if ev.inaxes == w.ax:
                    do_get = False
        if do_get:
            # resolve whether up happens in time-based plot
            # or a dynamic plot
            if ev.inaxes in self.dynamicPlots.values():
                self.getDynamicPoint(ev)
            else:
                self.getPoint(ev)
        self._mouseDrag = False


    def mouseMoveFn(self, ev):
        if self._mouseUp:
            self._mouseDrag = False
        else:
            self._mouseDrag = True


    def getPoint(self, ev):
        """
        If mouse click is released inside a time-based plot, change current time
        on time bar to that selected.
        """
        self.set_time(ev.xdata)


    def getDynamicPoint(self, ev):
        """
        If mouse clicks inside a user-specified 'dynamic' sub-plot axis,
        put the (x,y) coords of the click as a point onto the clipboard.
        """
        figName, layer_info, ixstr = self.subplot_lookup[ev.inaxes]
        fig_struct = self.plotter.figs[figName]
        # all layers share the same axes, so just get the first
        layer_struct = fig_struct.layers[layer_info[0]]
        self.clipboardPt = Point2D(ev.xdata, ev.ydata,
                                   xname=layer_struct.axes_vars[0],
                                   yname=layer_struct.axes_vars[1])
        print("Clipboard now contains: %s" % str(self.clipboardPt))


    def capturePoint(self, ev):
        """
        Create a dict of Point objects from the current time point in all sub-plots
        of all figures. Keys of the dict are the figure names
        Stored in the attribute capturedPts
        """
        pts_dict = {}
        # pts_dict : figName -> Point with coords t and <layer_names>
        # layer_dict : layer_name -> value

        for figName in self.plotter.figs:

            fig_struct = self.plotter.figs[figName]
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
                                pt_dict[data_name] = traj(self.t)['y']

            pts_dict.update({figName: Point(pt_dict)})

        if self.verbose >= 1:
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

    def goBack(self, ev):
        if self._last_ix is not None:
            self.updatePlots(self.times[self._last_ix])

    def updatePlots(self, new_time):
        self.set_time(new_time)
        self.plotter.updateDynamic(self.t, self.dynamicPlotFns)

    def refresh(self, ev):
        """For refresh button, e.g. use after zoom in dynamic plot
        """
        hard_reset = self._key_mod == 'shift'
        self.plotter.updateDynamic(self.t, self.dynamicPlotFns,
                                   hard_reset)

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


class context_object(object):
    # Abstract base class
    pass

class domain_GUI(context_object):
    pass # Not implemented yet!

class line_GUI(context_object):
    """
    Line of interest context_object for GUI
    """
    def __init__(self, gui, gui_axes, pt1, pt2):
        xnames = pt1.coordnames
        if pt2.coordnames != xnames:
            raise ValueError("Coordinate name mismatch")
        x1, y1 = pt1.coordarray
        x2, y2 = pt2.coordarray
        if x1 > x2:
            # ensure correct ordering for angles
            xt = x1
            yt = y1
            x1 = x2
            y1 = y2
            x2 = xt
            y2 = yt
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.dy = y2-y1
        self.dx = x2-x1
        self.length = np.linalg.norm((self.dx, self.dy))
        # angle relative to horizontal, in radians
        self.ang = atan2(self.dy,self.dx)
        self.ang_deg = 180*self.ang/pi
        # hook back to linked axes object in GUI
        self.gui_axes = gui_axes
        # declare self to GUI
        gui.declare_in_context(self)
        # move self to the currently selected object in GUI
        gui.selected_object = self
        self.gui = gui
        self.extra_fnspecs = {}
        self.extra_pars = {}
        self.extra_auxvars = {}
        self.extra_events = []
        print("Created line and moved to currently selected object")

        # actual MPL line object handle
        self.l = None
        self.name = '<untitled>'
        self.show()

    def __repr__(self):
        if self.extra_events == []:
            ev_str = '(no event)'
        else:
            ev_str = '(with event)'
        return "line_GUI(%.3f, %.3f, %.3f, %.3f) - '%s' %s" %(self.x1, self.y1, \
                                          self.x2, self.y2, self.name, ev_str)

    def show(self):
        if self.l is None:
            self.l = self.gui_axes.plot([self.x1,self.x2],
                                       [self.y1,self.y2],
                                   'y-')[0]
        else:
            self.l.set_visible(True)
        plt.draw()

    def unshow(self):
        self.l.set_visible(False)
        plt.draw()

    def remove(self):
        self.l.remove()
        plt.draw()

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

    def make_event_def(self, uniquename, dircode=0):
        self.name = uniquename
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
                    except Exception, e:
                        print("No attribute: '%s' in object in workspace '%s'" % (tracked_attr, wspace._name))
                        raise
                    plt.text(0.05, 0.05+i*0.04, '%s: %s = %.4g' % (obj_name, tracked_attr, data))
            plt.title('%s measures of %s (workspace: %s)'%(self.calc_context.sim.name, tracked_attr,
                                                           _escape_underscore(self.calc_context.workspace._name)))
            fig.canvas.set_window_title("Fig %i, Workspace %s" % (fignum, self.calc_context.workspace._name))
        plt.show()



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
                except Exception, e:
                    print("Failed to evaluate: '%s' in workspace '%s'" % (tracked.xstr, wspace._name))
                    raise
                try:
                    ydata = getattr(self.calc_context.workspace, tracked.ystr)
                except Exception, e:
                    print("Failed to evaluate: '%s' in workspace '%s'" % (tracked.ystr, wspace._name))
                    raise
                plt.plot(xdata, ydata,
                     tracked.style, label=_escape_underscore(tracked.ystr))
            plt.legend()
            plt.title('%s measures vs %s (workspace: %s)'%(self.calc_context.sim.name, tracked.xstr,
                                                           _escape_underscore(self.calc_context.workspace._name)))
            fig.canvas.set_window_title("Fig %i, Workspace %s" % (fignum, self.calc_context.workspace._name))
        plt.show()


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

global gui, tracker

# singleton pattern

# plotter is not a globally visible object outside of this module
plotter = plotter2D()

gui = diagnosticGUI(plotter)

tracker = tracker_manager()


