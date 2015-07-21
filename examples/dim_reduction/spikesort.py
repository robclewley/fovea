#from __future__ import division
from PyDSTool import *
import fovea
import fovea.domain2D as dom
from fovea import common, prep, graphics
#from fovea.graphics import gui, plotter
import ButterworthBandpass as bwb
import PyDSTool as dst
import PyDSTool.Toolbox.phaseplane as pp

from scipy.signal import butter, lfilter
import matplotlib as mpl
import matplotlib.pyplot as plt

class spikesorter(graphics.diagnosticGUI):
    def __init__(self, title):

        global plotter
        plotter = graphics.plotter2D()
        super().__init__(plotter)

        #Recover data:
        data = importPointset('shortsimdata1.dat',t=0,sep=',')

        vs = data['vararray'][0]
        vs = self.bandpass_filter(vs, 300, 3000, 32000)
        ts = data['t']

        self.traj = numeric_to_traj([vs], 'test_traj', ['x'], ts, discrete=False)

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

        self.setup({'11':
                   {'name': 'spikesort',
                    'scale': DOI,
                    'layers':['spikes'],
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

ssort = spikesorter("SSort")

cutoff = 15

ltarget = fovea.graphics.line_GUI(pp.Point2D(0, cutoff),
                                  pp.Point2D(15000, cutoff), subplot = '11')
ltarget.update(name ='threshline')

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

for num in result:
    plt.plot(num[0], cutoff, 'r*')

halt = True