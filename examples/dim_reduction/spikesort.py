#from __future__ import division
from PyDSTool import *
import fovea
import fovea.domain2D as dom
from fovea import common, prep, graphics
#from fovea.graphics import gui, plotter
import PyDSTool as dst
import PyDSTool.Toolbox.phaseplane as pp

from scipy.signal import butter, lfilter
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

        self.setup({'11':
                   {'name': 'spikesort',
                    'scale': DOI,
                    'layers':['spikes', 'thresh_crosses'],
                    'callbacks':'*',
                    'axes_vars': ['x', 'y']
                    }
                   },
                  size=(9, 7), with_times=False, basic_widgets=False)

        #Bad code carried over from fovea_game:
        fig_struct, figure = plotter._resolveFig(None)
        self.ax = fig_struct.arrange['11']['axes_obj']

        coorddict = {'x':
                     {'x':'t', 'layer':'spikes', 'style':'b.-'}
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

    def user_nav_func(self):
        cutoff =  ltarget.y1

        SS_event_args = {'name': 'SS_zerothresh',
                         'eventtol': 1e-3,
                         'eventdelay': 1e-4,
                         'starttime': 0,
                         'precise': True,
                         'active': True}
        #dircode = 1 is crossing from below
        SS_thresh_ev = Events.makePythonStateZeroCrossEvent('v', cutoff, 1, SS_event_args,
                                                            ssort.traj.variables['x'])

        #If I search the entire interval at once, it returns a partial list. Why?
        ts = ssort.traj.sample()['t']
        dt = ts[1]-ts[0]  # assumes uniformly timed samples
        result = SS_thresh_ev.searchForEvents((0, 15000), dt=dt, eventdelay=False)

        crosses = [num[0] for num in result]

        self.addDataPoints([crosses, [cutoff]*len(crosses)], layer='thresh_crosses', style='r*',
                           name='crossovers', force=True)
        plotter.show()

ssort = spikesorter("SSort")

cutoff = 20

ltarget = fovea.graphics.line_GUI(pp.Point2D(0, cutoff),
                                  pp.Point2D(15000, cutoff), subplot = '11')
ltarget.update(name ='threshline')

# For testing:

##class keyev(object):
##    pass
##
##ev = keyev
##ev.key = 'down'
##ssort.key_on(ev)

halt = True
