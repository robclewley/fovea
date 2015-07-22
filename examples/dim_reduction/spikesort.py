"""
shortsimdata1.data recovered from <http://www2.le.ac.uk/departments/engineering/research/bioengineering/neuroengineering-lab/software>

Based on work described in:
Martinez, J., Pedreira, C., Ison, M. J., Quian Quiroga, R. (2009). Realistic simulation of extracellular recordings.
In J Neurosci Methods, 184(2):285-93, doi: 10.1016/j.jneumeth.2009.08.017.
"""

#from __future__ import division
from PyDSTool import *
import fovea
import fovea.domain2D as dom
from fovea import common, prep, graphics
#from fovea.graphics import gui, plotter
import PyDSTool as dst
import PyDSTool.Toolbox.phaseplane as pp

from scipy.signal import butter, lfilter, argrelextrema
import matplotlib as mpl
import matplotlib.pyplot as plt

class spikesorter(graphics.diagnosticGUI):
    def __init__(self, title):

        global plotter
        plotter = graphics.plotter2D()
        graphics.diagnosticGUI.__init__(self, plotter)

        #Recover data:
        data = importPointset('shortsimdata1.dat',t=0,sep=',')

        vs = data['vararray'][0]
        vs = self.bandpass_filter(vs, 300, 3000, 32000)
        ts = data['t']

        self.traj = numeric_to_traj([vs], 'test_traj', ['x'], ts, discrete=False)

        locminx = argrelextrema(vs, np.less)[0]

        self.locminx = [x for x in locminx if self.traj(x)['x'] < 0]
        self.locminy = self.traj(self.locminx)['x']

        self.fovea_setup()

    def fovea_setup(self):
        #Setup code
        DOI = [(0,15000),(-30,30)]
        plotter.clean() # in case rerun in same session
        plotter.addFig('master',
                       title='spikesort',
                       xlabel='time', ylabel='mV',
                       domain=DOI)

        #Setup all layers
        plotter.addLayer('spikes')
        plotter.addLayer('thresh_crosses')
        plotter.addLayer('local_mins')
        plotter.addLayer('onsets')

        self.setup({'11':
                   {'name': 'spikesort',
                    'scale': DOI,
                    'layers':['spikes', 'thresh_crosses', 'local_mins', 'onsets'],
                    'callbacks':'*',
                    'axes_vars': ['x', 'y']
                    }
                   },
                  size=(9, 7), with_times=False, basic_widgets=False)

        #Bad code carried over from fovea_game:
        fig_struct, figure = plotter._resolveFig(None)
        self.ax = fig_struct.arrange['11']['axes_obj']

        coorddict = {'x':
                     {'x':'t', 'layer':'spikes', 'style':'b-'}
                     }
        self.addDataPoints(self.traj.sample(), coorddict = coorddict)
        self.addDataPoints([self.locminx, self.locminy], layer='local_mins', style='g*', name='mins', force= True)

        self.find_adjacent_mins(13793, 18.8)

        plotter.show()

    def model_namer(self):
        name = 'spikesorter_traj'
        return name

    def bandpass_filter(self, data, lowcut, highcut, fs, order= 5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype= 'band')
        y = lfilter(b, a, data)
        return y

    def user_nav_func(self):
        cutoff =  ltarget.y1

        SS_event_args = {'name': 'SS_zerothresh',
                         'eventtol': 1e-2,
                         'eventdelay': 1e-3,
                         'starttime': 0,
                         'active': True}
        #dircode = 1 is crossing from below
        SS_thresh_ev = Events.makePythonStateZeroCrossEvent('v', cutoff, 1, SS_event_args,
                                                            ssort.traj.variables['x'])

        #If I search the entire interval at once, it returns a partial list. Why?
        result = SS_thresh_ev.searchForEvents((0, 5000))
        result += SS_thresh_ev.searchForEvents((5000, 10000))
        result += SS_thresh_ev.searchForEvents((10000, 15000))

        crosses = [num[0] for num in result]

        self.addDataPoints([crosses, [cutoff]*len(crosses)], layer='thresh_crosses', style='r*', name='crossovers', force= True)
        plotter.show()

    def find_adjacent_mins(self, x, y):
        for i in range(0, len(self.locminx)):
            leftminx = self.locminx[i]
            rightminx = self.locminx[i + 1]

            leftminy = self.locminy[i]
            rightminy = self.locminy[i + 1]

            if x > leftminx and x < rightminx:
                break

        self.addDataPoints([leftminx, leftminy], style='y*', layer= 'onsets')
        self.addDataPoints([rightminx, rightminy], style='y*', layer= 'onsets')
        plotter.show()


ssort = spikesorter("SSort")

cutoff = 20

ltarget = fovea.graphics.line_GUI(ssort, pp.Point2D(0, cutoff),
                                  pp.Point2D(15000, cutoff), subplot = '11')
ltarget.update(name ='threshline')

class snap_point():
    def __init__(self, x, y):
        self.xdata = x
        self.ydata = y

sp = snap_point(4007, 17.26)
ssort.mouse_event_snap(sp)

halt = True