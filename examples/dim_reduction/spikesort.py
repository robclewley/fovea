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

#from neuro_data import * #CHANGE THIS IMPORT STATEMENT
from PyDSTool.Toolbox.neuro_data import *
from PyDSTool.Toolbox import data_analysis as da

from scipy.signal import butter, lfilter, argrelextrema
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

global default_sw
default_sw = 30

class spikesorter(graphics.diagnosticGUI):
    def __init__(self, title):

        #global plotter
        plotter = graphics.plotter2D()
        graphics.diagnosticGUI.__init__(self, plotter)

        #Recover data:
        data = importPointset('shortsimdata1.dat',t=0,sep=',')
        #data = importPointset('simdata1_100000.dat',t=0,sep=',')

        vs = data['vararray'][0]
        vs = self.bandpass_filter(vs, 300, 3000, 32000)
        ts = data['t']

        self.N = len(vs)
        self.traj = numeric_to_traj([vs], 'test_traj', ['x'], ts, discrete=False)

        self.selected_pcs = []
        self.fovea_setup()

    def fovea_setup(self):
        #Setup code
        DOI = [(0,self.N),(-30,30)]
        self.plotter.clean() # in case rerun in same session
        self.plotter.addFig('master',
                       title='spikesort',
                       xlabel='time', ylabel='mV',
                       domain=DOI)

        #Setup all layers
        self.plotter.addLayer('spikes')
        self.plotter.addLayer('thresh_crosses')
        self.plotter.addLayer('detected')
        self.plotter.addLayer('pcs')
        self.plotter.addLayer('scores')
        self.plotter.addLayer('loading_text', kind= 'text')

        self.setup({'11':
                   {'name': 'waveform',
                    'scale': DOI,
                    'layers':['spikes', 'thresh_crosses'],
                    'callbacks':'*',
                    'axes_vars': ['x', 'y']
                    },
                   '12':
                   {'name': 'detected spikes',
                    'scale': [(0, default_sw), (-80, 80)],
                    'layers':['detected', 'loading_text'],
                    #'callbacks':'*',
                    'axes_vars': ['x', 'y']
                    },
                   '21':
                    {'name': 'Principal Components',
                     'scale': [(0, default_sw), (-0.5, 0.5)],
                     'layers':['pcs'],
                     'callbacks':'*',
                     'axes_vars': ['x', 'y']
                     },
                    '22':
                    {'name': 'Classified Spikes',
                     'scale': [(-100, 100), (-100, 100)],
                     'layers':['scores'],
                     'callbacks':'*',
                     'axes_vars': ['x', 'y']
                     }
                   },
                  size=(12, 7), with_times=False, basic_widgets=True)

        #self.plotter.setText('load_perc', Loading: %d\%'%n, 'loading')
        self.plotter.addText(0.1, 0.8, 'Loading...', use_axis_coords=True,
                        name='load_perc', layer='loading_text', style='k')
        self.plotter.setLayer('loading_text', display=False)

        #Bad code carried over from fovea_game:
        fig_struct, figure = self.plotter._resolveFig(None)
        #self.ax = fig_struct.arrange['11']['axes_obj']

        coorddict = {'x':
                     {'x':'t', 'layer':'spikes', 'style':'b-'}
                     }
        self.addDataPoints(self.traj.sample(), coorddict = coorddict)

        evKeyOn = self.fig.canvas.mpl_connect('key_press_event', self.ssort_key_on)

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

    def user_pick_func(self, ev):
        if not ev.is_con_obj:
            fig_struct, figure = self.plotter._resolveFig(None)
            layer_struct = self.plotter._resolveLayer(figure, 'pcs')

            for name, artist in fig_struct['layers']['pcs']['handles'].items():
                if artist is ev.artist:
                    self.proj_PCs.insert(0, name)
                    self.proj_PCs = self.proj_PCs[0:2]

            for name in fig_struct['layers']['pcs']['handles'].keys():
                if name not in self.proj_PCs:
                    #fig_struct['layers']['pcs']['handles'][name].set_linewidth(1)
                    #self.plotter.setData('pcs', name=name, linewidth= 1)
                    layer_struct.data[name].update({'linewidth': 1})

            for pc in self.proj_PCs:
                #self.plotter.setData('pcs', name= pc, linewidth= 2.5)
                #fig_struct['layers']['pcs']['handles'][pc].set_linewidth(2.5)
                layer_struct.data[pc].update({'linewidth': 2.5})

            self.plotter.show()

            self.proj_vec1 = fig_struct['layers']['pcs']['handles'][self.proj_PCs[0]].get_ydata()
            self.proj_vec2 = fig_struct['layers']['pcs']['handles'][self.proj_PCs[1]].get_ydata()

            self.project_to_PC()

            #ev.artist.set_linewidth(2.5)

            #for pc, pc_artist in fig_struct['layers']['pcs']['handles'].items():
                #if pc_artist is ev.artist:
                    #if not pc_artist.get_linewidth() == 2.5:
                        #odd_vec = self.proj_vec2
                        #self.proj_vec2 = self.proj_vec1
                        #self.proj_vec1 = fig_struct['layers']['pcs']['data'][pc]['data'][1]

                        #pc_artist.set_linewidth(2.5)
                #else:
                    #pc_artist.set_linewidth(1)

            #self.project_to_PC()


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
            #result = SS_thresh_ev.searchForEvents((0, 15000), dt=dt, eventdelay=False)
            result = SS_thresh_ev.searchForEvents((0, self.N), dt=dt, eventdelay=False)

            crosses = [num[0] for num in result]

            self.plotter.setLayer('loading_text', display= True)

            self.addDataPoints([crosses, [cutoff]*len(crosses)], layer='thresh_crosses', style='r*', name='crossovers', force= True)
            self.X = self.compute_bbox(crosses)

            #Plot detected spike
            print('self.X.shape', self.X.shape)
            if len(self.X.shape) == 1:
                self.addDataPoints([list(range(0, len(self.X))), self.X], layer= 'detected', style= 'b-', name= 'spike0', force= True)

            else:
                c= 0
                for spike in self.X:
                    self.addDataPoints([list(range(0, len(spike))), spike], layer= 'detected', style= 'b-', name= 'spike'+str(c), force= True)
                    c += 1

            self.plotter.show()

    def compute_bbox(self, crosses):
        try:
            search_width = self.context_objects['ref_box'].dx
            self.context_objects['ref_box'].remove()
        except KeyError:
            print("No 'ref_box' defined. Defaulting spike search width to",default_sw)
            search_width = default_sw

        fig_struct, figs = self.plotter._resolveFig(None)

        #Clear existing bounding boxes
        rem_names = []
        for con_name, con_obj in self.context_objects.items():
            if isinstance(con_obj, box_GUI) and con_name is not 'ref_box':
                rem_names.append(con_name)
        for name in rem_names:
            self.context_objects[name].remove(draw= False)
            #try:
                #del fig_struct['layers']['detected']['data']['det_'+name]
            #except KeyError:
                #pass

        fig_struct['layers']['detected']['data'] = {}

        self.plotter.show(rebuild= True)

        #Create new bounding boxes
        c = 0
        for x in crosses:
            thresh_ix = round(x)
            tlo = thresh_ix - round(search_width/2)
            thi = thresh_ix + round(search_width/2)
            result = find_internal_extrema(self.traj.sample(tlo= tlo, thi= thi))

            #Center box around spike peak
            tlo = tlo + result['global_max'][0] - math.floor(search_width/2)
            #thi = tlo + result['global_max'][0] + round(search_width/2)
            thi = tlo + search_width
            box_GUI(self, pp.Point2D(tlo, result['global_max'][1]), pp.Point2D(thi, result['global_min'][1]),
                        name= 'spike'+str(c), select= False)

            #Pin data
            coorddict = {'x':
                         {'x':'t', 'layer':'detected', 'style':'b-', 'name':'det_spike'+str(c)}
                         }
            ##ISSUE: spike_seg length becomes very small for any self.x1 > 100000.
            spike_seg = ssort.context_objects['spike'+str(c)].pin_contents(self.traj, coorddict)

            #print('length of spike_seg:',len(spike_seg))
            #print('c:',c)
            try:
                X = np.row_stack((X, spike_seg))
            except NameError:
                X = spike_seg
            #except ValueError: #Returns wrong concatenation size.
                #pass

            self.plotter.setText('load_perc', 'Complete: %0.2f'%(c/len(crosses)), 'loading_text')

            c += 1

        self.plotter.setLayer('loading_text', display=False)
        self.set_selected_object(self.context_objects['threshline'])

        return X

    def project_to_PC(self):
        print('shape proj vecs:', np.column_stack((self.proj_vec1, self.proj_vec2)).shape)
        print('shape X', self.X.shape)
        Y = np.dot(self.X, np.column_stack((self.proj_vec1, self.proj_vec2)))

        #Y = self.p._execute(X, 2)
        self.addDataPoints([Y[:,0], Y[:,1]], layer='scores', style='k*', name='lo_D', force= True)
        self.show(rebuild = True)

    def ssort_key_on(self, ev):
        self._key = k = ev.key  # keep record of last keypress
        fig_struct, fig = self.plotter._resolveFig(None)

        if k == 'p':
            print('doing PCA...')
            try:
                X = self.X
            except AttributeError:
                print('Must detect spikes before performing PCA.')
                return

            self.p = da.doPCA(X, len(X[0]), len(X[0]))
            self.proj_vec1 = self.p.get_projmatrix()[:, 0]
            self.proj_vec2 = self.p.get_projmatrix()[:, 1]

            print('proj mat shape pre-transpose: ', self.p.get_projmatrix().shape)
            #Y3 = self.p.get_projmatrix()
            #Y2 = p._execute(X, 2)
            #Y3 = p._execute(X, 4)

            self.addDataPoints([list(range(0, len(self.proj_vec1))) , self.proj_vec1], style= 'r-', layer= 'pcs', linewidth = 2.5, name= 'first_pc', force= True)
            self.addDataPoints([list(range(0, len(self.proj_vec2))) , self.proj_vec2], style= 'g-', layer= 'pcs', linewidth = 2.5, name= 'second_pc', force= True)

            self.plotter.show()

            ##ISSUE: This is a very inconvenient way to set a linewidth
            self.proj_PCs = ['first_pc', 'second_pc']

            try:
                self.addDataPoints([list(range(0, len(self.p.get_projmatrix()))) ,self.p.get_projmatrix()[:,2]],
                                   style= 'y-', layer= 'pcs', name= 'third_pc', force= True)
            except IndexError:
                pass

            self.show()

            self.project_to_PC()

ssort = spikesorter("SSort")

cutoff = 20

rets = find_internal_extrema(ssort.traj.sample())

fig_struct, figs = ssort.plotter._resolveFig(None)

halt = True

##class keyev(object):
##    pass
##
##ev = keyev
##ev.key = 'down'
##ssort.key_on(ev)

ssort.show_tree()

halt = True
