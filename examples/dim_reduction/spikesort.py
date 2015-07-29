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
from fovea.graphics import *
#from fovea.graphics import gui, plotter
import PyDSTool as dst
import PyDSTool.Toolbox.phaseplane as pp

from neuro_data import *

from scipy.signal import butter, lfilter, argrelextrema
import matplotlib as mpl
import matplotlib.pyplot as plt

class spikesorter(graphics.diagnosticGUI):
    def __init__(self, title):

        #global plotter
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
        self.plotter.clean() # in case rerun in same session
        self.plotter.addFig('master',
                       title='spikesort',
                       xlabel='time', ylabel='mV',
                       domain=DOI)

        #Setup all layers
        self.plotter.addLayer('spikes')
        self.plotter.addLayer('thresh_crosses')
        self.plotter.addLayer('detected')

        self.setup({'11':
                   {'name': 'waveform',
                    'scale': DOI,
                    'layers':['spikes', 'thresh_crosses'],
                    'callbacks':'*',
                    'axes_vars': ['x', 'y']
                    },
                   '12':
                   {'name': 'detected spikes',
                    'scale': [(0, 75), (-30, 30)],
                    'layers':['detected'],
                    #'callbacks':'*',
                    'axes_vars': ['x', 'y']
                    }
                   },
                  size=(9, 7), with_times=False, basic_widgets=False)

        #Bad code carried over from fovea_game:
        fig_struct, figure = self.plotter._resolveFig(None)
        #self.ax = fig_struct.arrange['11']['axes_obj']

        coorddict = {'x':
                     {'x':'t', 'layer':'spikes', 'style':'b.-'}
                     }
        self.addDataPoints(self.traj.sample(), coorddict = coorddict)

        self.plotter.show()

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
        if self.selected_object.name is 'threshline':
            cutoff =  self.selected_object.y1

            SS_event_args = {'name': 'SS_zerothresh',
                             'eventtol': 1e-3,
                             'eventdelay': 1e-4,
                             'starttime': 0,
                             'precise': True,
                             'active': True}
            #dircode = 1 is crossing from below
            SS_thresh_ev = Events.makePythonStateZeroCrossEvent('v', cutoff, 1, SS_event_args,
                                                                var= ssort.traj.variables['x'])

            #If I search the entire interval at once, it returns a partial list. Why?
            ts = ssort.traj.sample()['t']
            dt = ts[1]-ts[0]  # assumes uniformly timed samples
            result = SS_thresh_ev.searchForEvents((0, 15000), dt=dt, eventdelay=False)

            crosses = [num[0] for num in result]

            self.addDataPoints([crosses, [cutoff]*len(crosses)], layer='thresh_crosses', style='r*', name='crossovers', force= True)
            self.compute_bbox(crosses)

        #plotter.show()

    def compute_bbox(self, crosses):
        try:
            search_width = self.context_objects['ref_box'].dx
            self.context_objects['ref_box'].remove()
        except KeyError:
            print("No 'ref_box' defined. Defaulting spike search width to 75.")
            search_width = 75


        fig_struct, figs = self.plotter._resolveFig(None)

        #Clear existing bounding boxes
        rem_names = []
        for con_name, con_obj in self.context_objects.items():
            if isinstance(con_obj, box_GUI) and con_name is not 'ref_box':
                rem_names.append(con_name)
        for name in rem_names:
            self.context_objects[name].remove(draw= False)
            del fig_struct['layers']['detected']['data']['det_'+name]

        self.plotter.show(rebuild= True)

        #Create new bounding boxes
        c = 0
        for x in crosses:
            thresh_ix = round(x)
            tlo = thresh_ix - round(search_width/2)
            thi = thresh_ix + round(search_width/2)
            result = find_internal_extrema(self.traj.sample(tlo= tlo, thi= thi))

            #Center box around spike peak
            tlo = tlo + result['global_max'][0] - round(search_width/2)
            thi = tlo + result['global_max'][0] + round(search_width/2)
            box_GUI(self, pp.Point2D(tlo, result['global_max'][1]), pp.Point2D(thi, result['global_min'][1]),
                        name= 'spike'+str(c),select= False)

            #Pin data
            coorddict = {'x':
                         {'x':'t', 'layer':'detected', 'style':'b-', 'name':'det_spike'+str(c)}
                         }
            ssort.context_objects['spike'+str(c)].pin_contents(self.traj, coorddict)
            c += 1

        self.set_selected_object(self.context_objects['threshline'])

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

        self.addDataPoints([crosses, [cutoff]*len(crosses)], layer='thresh_crosses', style='r*',
                           name='crossovers', force=True)
        plotter.show()

ssort = spikesorter("SSort")

cutoff = 20
#ltarget = fovea.graphics.line_GUI(ssort, pp.Point2D(0, cutoff),
                                  #pp.Point2D(15000, cutoff), subplot = '11')
#ltarget.update(name ='threshline')

rets = find_internal_extrema(ssort.traj.sample())

fig_struct, figs = ssort.plotter._resolveFig(None)

halt = True

#class snap_point():
    #def __init__(self, x, y):
        #self.xdata = x
        #self.ydata = y

#sp = snap_point(-0.51, 0.51)
#ssort.mouse_event_snap(sp)
# For testing:

##class keyev(object):
##    pass
##
##ev = keyev
##ev.key = 'down'
##ssort.key_on(ev)

halt = True
