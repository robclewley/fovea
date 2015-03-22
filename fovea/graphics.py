"""
Data visualization class for dynamical systems

Bryce Chung and Rob Clewley, 2012


==== Usage notes:
Plotting styles can be given either as a string or as a dictionary of
  style kwargs suitable for plot command

"""


import matplotlib.pyplot as pp
from matplotlib.widgets import Slider, Button
import numpy as np
from copy import copy

from PyDSTool import args, numeric_to_traj, Point
from PyDSTool.Toolbox.phaseplane import Point2D


# ----------------------------------------------

class plotter2D(object):

    colors = ['b', 'g', 'r', 'c', 'm', 'k', 'y']

    def __init__(self):
        self.figs = {}
        self._max_fig_num = 0
        self.currFig = None

    def clean(self):
        """
        Delete all figure data (doesn't clear figures)
        """
        self.figs = {}
        self._max_fig_num = 0

    def autoScale(self, figure=None):
        xScale = [0,1]
        yScale = [0,1]

        for figName, fig in self.figs.iteritems():
            for layerName, layer in fig.layers.iteritems():
                for dName, d in layer['data'].iteritems():
                    xScale[0] = min(min(d[0][0]), xScale[0])
                    xScale[1] = max(max(d[0][0]), xScale[1])
                    yScale[0] = min(min(d[0][1]), yScale[0])
                    yScale[1] = max(max(d[0][1]), yScale[1])

        pp.xlim(xScale)
        pp.ylim(yScale)


    # This method is never called!
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
##            fig = pp.figure(fig_struct.fignum)
##            #pp.title(fig_struct.title)
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
##                for ix, subplot in fig_struct.arrange.iteritems():
##                    ax = fig.add_subplot(shape[0], shape[1], shape[1]*(int(ix[0])-1)+int(ix[1]))
##                    try:
##                        axesLabels = subplot['axes']
##                        if axesLabels is not None:
##                            ax.set_xlabel(axesLabels[0])
##                            ax.set_ylabel(axesLabels[1])
##                    except:
##                        raise ValueError("Error setting axis label for subplot")
##                    try:
##                        pp.title(subplot['name'])
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
##                    pp.autoscale(enable=autoscale)
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

        for ixstr, spec in arrPlots.iteritems():
            if subplot is not None and subplot not in \
                              [ixstr, spec['name']]:
                continue
            layer_info = spec['layers']
            for layer in layer_info:
                lay = fig_struct.layers[layer]
                if not lay.display or lay.kind != 'data':
                    continue
                print "\nLAYER:", layer
                print "  sub-plot: ", ixstr
                print "  style", lay.style
                print "  axes:", lay.axes_vars
                print "  data:"
                for dname, dstruct in lay.data.items():
                    if not dstruct['display']:
                        continue
                    print "     name:", dname, " style:", dstruct['style']



    ## Figure Management Tools ##

    def addFig(self, label, title="", xlabel="", ylabel="", tdom=None, display=True):
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

        for ixstr, spec in arrPlots.iteritems():
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

        # Check to see layer does not already exist #
        if fig_struct.layers.has_key(layer_name):
            raise KeyError("Layer name already exists in figure!")

        # default plot style generated from list of colors
        color = self.colors[len(fig_struct.layers)%len(self.colors)]
        line = '-'
        style = color+line

        layAttrs = args()
        layAttrs.data = {}
        layAttrs.display = True
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
            # Check to see that parameter exists in layers #
            ## Possibly change to account for different properties of specific artist objects?
            if not layAttrs.has_key(kw):
                raise KeyError("Parameter is not a property of the layer.")

            layAttrs[kw] = kwargs[kw]

        fig_struct.layers[layer_name] = layAttrs


    def setLayer(self, label, figure=None, **kwargs):
        """
        Allows users to arrange data sets in a figure
        """
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
            raise KeyError("Data name already exists in layer, %s" %name)

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

        #print "Adding: %s >> %s" % (figure, name)
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
                        #print "updateDynamic calling function:", dynamicFns[layer]
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

        for dkey, dstruct in lay.data.iteritems():
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


def castNull(null):
    """
    Utility to convert nullcline object to a set of points useable by plotter.
    """
    return [null.array[:,0], null.array[:,1]]



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


    def buildPlotter2D(self, figsize=None):
        """
        Create time bar widget.
        Create capture points widget.
        Also closes all figure windows and resets masterWin attribute.

        Optional size is a pair of figure screen size measurements in inches
        """

        pp.close('all')
        self.masterWin = None

        for figName, fig_struct in self.plotter.figs.items():
            if figsize is not None:
                fig_handle = pp.figure(fig_struct.fignum, figsize=figsize)
                # fig sizing later, with fig.set_figsize_inches doesn't seem
                # to work properly unless put in the initialization call
            else:
                fig_handle = pp.figure(fig_struct.fignum)
            fig_handle.canvas.set_window_title(fig_struct.title + " : Master window")
            if figName != 'Master':
                continue
            ##### Set up master window controls
            pp.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.1,
                               wspace=0.2, hspace=0.23)

            self.masterWin = fig_handle

            ## Time bar controls time lines in figures
            sliderRange = self.times
            slide = pp.axes([0.25, 0.02, 0.65, 0.03])
            tMin = min(sliderRange)
            tMax = max(sliderRange)
            if self.t is None:
                self.set_time( (tMin + tMax)/2. )
            self.widgets['timeBar'] = Slider(slide, 'Time', tMin, tMax,
                                            valinit=self.t, color='r',
                                            dragging=False, valfmt='%1.4f')

            # button axes are in figure coords: (left, bottom, width, height)

            ## +/- dt buttons
            m_dt_Button = Button(pp.axes([0.16, 0.02, 0.017, 0.03]), '-dt')
            self.widgets['minus_dt'] = m_dt_Button

            p_dt_Button = Button(pp.axes([0.18, 0.02, 0.017, 0.03]), '+dt')
            self.widgets['plus_dt'] = p_dt_Button

            ## Capture point button in lower left
            captureButton = Button(pp.axes([0.055, 0.02, 0.08, 0.03]), 'Capture Point')
            self.widgets['capturePoint'] = captureButton

            ## Refresh button
            refreshButton = Button(pp.axes([0.005, 0.02, 0.045, 0.03]), 'Refresh')
            self.widgets['refresh'] = refreshButton

            ## Go back to last point button
            backButton = Button(pp.axes([0.005, 0.06, 0.045, 0.03]), 'Back')
            self.widgets['goBack'] = backButton

            ## Build up each subplot, left to right, top to bottom
            shape = fig_struct.shape
            for i in range(shape[0]):
                for j in range(shape[1]):
                    ixstr = str(i+1) + str(j+1)
                    try:
                        subplot_struct = fig_struct.arrange[ixstr]
                    except KeyError:
                        continue
                    layer_info = subplot_struct['layers']
                    if not isinstance(layer_info, list):
                        # singleton string layer name
                        layer_info = [layer_info]

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

            fig_handle.canvas.draw()

        # activate callbacks
        self.widgets['timeBar'].on_changed(self.updatePlots)
        self.widgets['capturePoint'].on_clicked(self.capturePoint)
        self.widgets['refresh'].on_clicked(self.refresh)
        self.widgets['minus_dt'].on_clicked(self.minus_dt)
        self.widgets['plus_dt'].on_clicked(self.plus_dt)
        self.widgets['goBack'].on_clicked(self.goBack)

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
            for ixstr, spec in arrPlots.iteritems():
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
        print "Clipboard now contains:", self.clipboardPt


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
                print "\n\n=============================="
                print "Figure: %s" % figName
                print "@ time = %.3f" % pt['t']
                print pt
            print "\n"

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
            pp.draw()

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


# -------------------

def checked_scale(sc):
    """Internal utility to verify scale syntax:
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


# ---------------------------------------------------------

global gui

# plotter is not a globally visible object outside of this module
plotter = plotter2D()

gui = diagnosticGUI(plotter)

