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

from PyDSTool.Toolbox.neuro_data import *
#from PyDSTool.Toolbox import data_analysis as da
from mdp.nodes import PCANode

from scipy.signal import butter, lfilter, argrelextrema
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import time

global default_sw
default_sw = 64

global show_after_load
show_after_load = True

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

        #Threshold used in martinez et al.
        self.mthresh = 4*(np.median(np.abs(vs))/0.6475)

        self.selected_pcs = []
        self.fovea_setup()

        #x = self.traj.sample()['x']
        #r = x > 10
        #above_thresh = np.where(r == True)[0]

        #spike = []
        #peaks = []
        #crosses = []
        #last_i = above_thresh[0] - 1
        #for i in above_thresh:
            #if i - 1 != last_i:
                #crosses.append(i)
                ##Return x value of the highest y value.
                #peaks.append(np.where(x == max(list(x[spike])))[0][0])
                #spike = []
            #spike.append(i)
            #last_i = i

        print("STEP 1:")
        print("Create a horizontal line of interest by pressing 'l'.")
        print("Once created, this line can be forced to extent by pressing 'm'.")
        print("Enter 'ssort.selected_object.update(name = 'thresh')' to identify the line as a threshold for spike detection")
        print("Once the line is renamed to 'thresh', the arrow keys can be used to move it up and down.")
        self.tutorial = 'step2'

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
                     #'scale': [(-100, 100), (-100, 100)],
                     'scale': [(-300, 300), (-300, 300)],
                     'layers':['scores'],
                     'callbacks':'*',
                     'axes_vars': ['x', 'y']
                     }
                   },
                  size=(8, 8), with_times=False, basic_widgets=True)

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

        self.plotter.auto_scale_domain(subplot= '11', xcushion= 0)

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

            #Clicked a PC.
            if ev.mouseevent.inaxes is fig_struct['arrange']['21']['axes_obj']:
                layer_struct = self.plotter._resolveLayer(figure, 'pcs')

                for name, artist in fig_struct['layers']['pcs']['handles'].items():
                    if artist is ev.artist:
                        self.proj_PCs.insert(0, name)
                        self.proj_PCs = self.proj_PCs[0:2]

                for name in fig_struct['layers']['pcs']['handles'].keys():
                    if name not in self.proj_PCs:
                        #fig_struct['layers']['pcs']['handles'][name].set_linewidth(1)
                        #self.plotter.setData('pcs', name=name, linewidth= 1)
                        ##ISSUE: This a very unintuitive way to update width, but it wont stick otherwise.
                        #layer_struct.data[name].update({'linewidth': 1})
                        self.plotter.setData2(name, layer='pcs', linewidth= 1)

                for pc in self.proj_PCs:
                    #self.plotter.setData('pcs', name= pc, linewidth= 2.5)
                    #fig_struct['layers']['pcs']['handles'][pc].set_linewidth(2.5)
                    #layer_struct.data[pc].update({'linewidth': 2.5})
                    self.plotter.setData2(pc, layer='pcs', linewidth= 2.5)

                self.plotter.show()

                self.proj_vec1 = fig_struct['layers']['pcs']['handles'][self.proj_PCs[0]].get_ydata()
                self.proj_vec2 = fig_struct['layers']['pcs']['handles'][self.proj_PCs[1]].get_ydata()

                self.project_to_PC()

            #Clicked a detected spike
            picked_spikes = None
            if ev.mouseevent.inaxes is fig_struct['arrange']['12']['axes_obj']:
                picked_spikes = 'detected'
            if ev.mouseevent.inaxes is fig_struct['arrange']['22']['axes_obj']:
                picked_spikes = 'scores'

            if picked_spikes is not None:
                detected_struct = self.plotter._resolveLayer(figure, 'detected')
                scores_struct = self.plotter._resolveLayer(figure, 'scores')

                for name, artist in fig_struct['layers'][picked_spikes]['handles'].items():
                    if artist is ev.artist:
                        #Saving these variables for convenience.
                        self.pick_det = {name : detected_struct.data[name]}

                        #detected_struct.data[name].update({'linewidth': 3, 'zorder':10, 'style':'y-'})
                        self.plotter.setData2(name, layer='detected', linewidth= 3, zorder= 10, style='y-')
                        try:
                            self.pick_score = scores_struct.data[name]
                            self.plotter.setData2(name, layer='scores', markersize= 12, zorder= 10, style='y*')
                            #scores_struct.data[name].update({'linewidth': 3, 'zorder':10,
                                                             #'markersize': 12, 'style':'y*'})
                        except KeyError:
                            pass

                    else:
                        #detected_struct.data[name].update({'linewidth': 1, 'zorder':1, 'style':self.default_colors[name]+str('-')})
                        self.plotter.setData2(name, layer='detected', linewidth= 1, zorder= 1, style= self.default_colors[name]+str('-'))
                        try:
                            self.plotter.setData2(name, layer='scores', markersize= 6, zorder= 1, style= self.default_colors[name]+str('*'))
                            #scores_struct.data[name].update({'linewidth': 1, 'zorder':1,
                                                             #'markersize': 6, 'style':self.default_colors[name]+str('*')})
                        except KeyError:
                            pass

                ##ISSUE: Have to rebuild in order to update data attributes.
                self.plotter.show(rebuild= True)


    def user_nav_func(self):
        if self.selected_object.name is 'thresh':
            try:
                self.search_width = self.context_objects['ref_box'].dx
                self.context_objects['ref_box'].remove()
            except KeyError:
                self.search_width = default_sw

            if self.tutorial == 'step2':
                print("STEP 2: ")
                print("When thresh is in place, press 'd' to capture each spike crossing the threshold in a bounding box.")
                print("Each detected spike will be placed in the top right subplot.")
                self.tutorial = 'step3'

            self.plotter.setText('load_perc', 'Loading...', 'loading_text')
            self.plotter.setLayer('loading_text', display= True)
            self.show()

            cutoff =  self.selected_object.y1

            #SS_event_args = {'name': 'SS_zerothresh',
                             #'eventtol': 1e-3,
                             #'eventdelay': 1e-4,
                             #'eventinterval': self.search_width,
                             #'starttime': 0,
                             #'precise': True,
                             #'active': True}
            ##dircode = 1 is crossing from below
            #SS_thresh_ev = Events.makePythonStateZeroCrossEvent('v', cutoff, 1, SS_event_args,
                                                                #var= ssort.traj.variables['x'])

            ##If I search the entire interval at once, it returns a partial list. Why?
            #ts = ssort.traj.sample()['t']
            #dt = ts[1]-ts[0]  # assumes uniformly timed samples
            ##dt = default_sw
            ##result = SS_thresh_ev.searchForEvents((0, 15000), dt=dt, eventdelay=False)
            #result = SS_thresh_ev.searchForEvents((0, self.N), dt=dt, eventdelay=False)
            #r = ssort.traj.sample()['x'] > cutoff

            x = self.traj.sample()['x']
            r = x > cutoff
            above_thresh = np.where(r == True)[0]

            spike = []
            peaks = []
            crosses = []

            last_i = above_thresh[0] - 1
            for i in above_thresh:
                if i - 1 != last_i:
                    crosses.append(i)
                    #Return x value of the highest y value.
                    peaks.append(np.where(x == max(list(x[spike])))[0][0])
                    spike = []
                spike.append(i)
                last_i = i

            self.crosses = crosses
            self.peaks = peaks
            #self.spikes = spikes

            #self.addDataPoints([crosses, [cutoff]*len(crosses)], layer='thresh_crosses', style='r*', name='crossovers', force= True)
            self.plotter.addData([self.crosses, [cutoff]*len(self.crosses)], layer='thresh_crosses', style='r*', name='crossovers', force= True)

            self.plotter.setLayer('loading_text', display= False)
            self.show()


    def compute_bbox(self):

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
        for i in range(len(self.peaks)):
            #thresh_ix = round(x)
            #tlo = thresh_ix - round(self.search_width/2)
            #thi = thresh_ix + round(self.search_width/2)
            #result = find_internal_extrema(self.traj.sample(tlo= tlo, thi= thi))

            ##Center box around spike peak
            ##tlo = tlo + result['global_max'][0] - math.floor(search_width/2)

            ##Place peak at step 20 (as in Martinez et al.)
            #tlo = tlo + result['global_max'][0] - 20
            #thi = tlo + self.search_width
            tlo = self.peaks[i] - 20
            thi = tlo + self.search_width
            valley = min(self.traj.sample()['x'][tlo:thi])
            box_GUI(self, pp.Point2D(tlo, self.traj.sample()['x'][self.peaks[i]]),
                    pp.Point2D(thi, valley),name= 'spike_box'+str(c), select= False)


            #Pin data
            coorddict = {'x':
                         {'x':'t', 'layer':'detected', 'style':'k-', 'name':'det_spike'+str(c)}
                         }
            ##ISSUE: spike_seg length becomes very small for any self.x1 > 100000.
            spike_seg = ssort.context_objects['spike_box'+str(c)].pin_contents(self.traj, coorddict)

            print('length of spike_seg:',len(spike_seg))
            print('c:',c)
            try:
                X = np.row_stack((X, spike_seg))
            except NameError:
                X = spike_seg
            #except ValueError: #Returns wrong concatenation size.
                #pass

            self.plotter.setText('load_perc', 'Complete: %0.2f'%(c/len(self.crosses)), 'loading_text')
            if show_after_load:
                self.show()

            c += 1

        self.plotter.setLayer('loading_text', display=False)
        self.set_selected_object(self.context_objects['thresh'])

        return X

    def project_to_PC(self):
        #print('shape proj vecs:', np.column_stack((self.proj_vec1, self.proj_vec2)).shape)
        #print('shape X', self.X.shape)
        Y = np.dot(self.X, np.column_stack((self.proj_vec1, self.proj_vec2)))

        #If moving to a smaller number of spikes, just forcing out data by reassigning names won't work. Must clear.
        self.clearData('scores')
        self.show()

        self.default_colors = {}
        #Add spikes as individual lines, so they can be referenced individually.
        c = 0
        for spike in Y:
            name = 'spike'+str(c)
            self.default_colors[name] = 'k'
            self.addDataPoints([spike[0], spike[1]], layer='scores', style=self.default_colors[name]+'*', name= name)
            c += 1

        self.plotter.auto_scale_domain(subplot = '22')

        self.show(rebuild = True)

    def ssort_key_on(self, ev):
        self._key = k = ev.key  # keep record of last keypress
        fig_struct, fig = self.plotter._resolveFig(None)

        class_keys = ['1','2','3','4']

        if k in class_keys:
            if isinstance(self.selected_object, box_GUI):
                for dname, dstruct in fig_struct['layers']['scores']['data'].items():
                    if self.selected_object.x1 < dstruct['data'][0] < self.selected_object.x2 and \
                    self.selected_object.y1 < dstruct['data'][1] < self.selected_object.y2:
                        if k == '1':
                            self.default_colors[dname] = 'r'
                            self.plotter.setData2(dname, layer='detected', style= 'r-')
                            self.plotter.setData2(dname, layer='scores', style= 'r*')
                        if k == '2':
                            self.default_colors[dname] = 'g'
                            self.plotter.setData2(dname, layer='detected', style= 'g-')
                            self.plotter.setData2(dname, layer='scores', style= 'g*')

                        if k == '3':
                            self.default_colors[dname] = 'b'
                            self.plotter.setData2(dname, layer='detected', style= 'b-')
                            self.plotter.setData2(dname, layer='scores', style= 'b*')

            #ISSUE: This only works when rebuild true.
            self.plotter.show(rebuild = True)

        if k== 'd':
            try:
                self.crosses
            except AttributeError:
                print("Can't detect spikes until threshold crossings have been found.")
                return

            self.plotter.setLayer('loading_text', display=True)
            self.X = self.compute_bbox()

            self.default_colors = {}

            if len(self.X.shape) == 1:
                self.default_colors['spike0'] = 'k'
                self.addDataPoints([list(range(0, len(self.X))), self.X], layer= 'detected', style= self.default_colors['spike0']+'-', name= 'spike0', force= True)

            else:
                c= 0
                for spike in self.X:
                    name = 'spike'+str(c)
                    self.default_colors[name] = 'k'
                    self.addDataPoints([list(range(0, len(spike))), spike], layer= 'detected', style= self.default_colors[name]+'-', name= name, force= True)
                    c += 1

            self.plotter.auto_scale_domain(xcushion = 0, subplot = '12')
            self.show()

            if self.tutorial == 'step3':
                print("STEP 3: ")
                print("You can now press 'p' to perform PCA on the detected spikes.")
                print("The bottom right subplot will display the first 3 principal components (in red, green, and yellow respectively.)")
                print("The bottom left subplot will show the detected spikes projected onto the first two PCs")
                self.tutorial = 'step4'

        if k == 'p':
            print('doing PCA...')
            try:
                X = self.X
            except AttributeError:
                print('Must detect spikes before performing PCA.')
                return

            #self.p = da.doPCA(X, len(X[0]), len(X[0]))
            self.p = PCANode(output_dim=0.99, reduce= True, svd= True)
            self.p.train(X)
            self.proj_vec1 = self.p.get_projmatrix()[:, 0]
            self.proj_vec2 = self.p.get_projmatrix()[:, 1]

            #print('proj mat shape: ', self.p.get_projmatrix().shape)

            self.addDataPoints([list(range(0, len(self.proj_vec1))) , self.proj_vec1], style= 'r-', layer= 'pcs', linewidth = 2.5, name= 'first_pc', force= True)
            self.addDataPoints([list(range(0, len(self.proj_vec2))) , self.proj_vec2], style= 'g-', layer= 'pcs', linewidth = 2.5, name= 'second_pc', force= True)

            self.plotter.show()
            self.proj_PCs = ['first_pc', 'second_pc']

            try:
                self.addDataPoints([list(range(0, len(self.p.get_projmatrix()))) ,self.p.get_projmatrix()[:,2]],
                                   style= 'y-', layer= 'pcs', name= 'third_pc', force= True)
            except IndexError:
                pass

            self.plotter.auto_scale_domain(xcushion = 0, subplot = '21')
            self.show()

            self.project_to_PC()

            if self.tutorial == 'step4':
                print("STEP 4: ")
                print("Use mouse clicks to explore the data.")
                print("Clicking on detected spikes in the top-right will highlight the corresponding projection in the bottom right (and vice versa).")
                print("You can also change the set of PCs onto which the data are projected by clicking the desired project PCs in the bottom left")
                print("Try exploring simdata1_100000.dat, which is a larger dataset and reveals clearer spiking clusters.")

                print("NOTE ALSO: ")
                print("Creating a bounding box in the upper-left plot and renaming it to 'ref_box', will change the search width of the detected spike.")
                print("e.g., if you want detected spikes to be 30 msec long, the box's .dx value must be 30.")
                print("After creating the box, it will be set to the current selected object. You can select the thresh line again by clicking on it.")


ssort = spikesorter("SSort")

rets = find_internal_extrema(ssort.traj.sample())

fig_struct, figs = ssort.plotter._resolveFig(None)

halt = True

##threshline = line_GUI(ssort, pt1, pt2)
#threshline = line_GUI(ssort, pp.Point2D(0,ssort.mthresh), pp.Point2D(100000, ssort.mthresh), layer='gx_objects', subplot='11',
                     #name='thresh', select=True)
#ssort.navigate_selected_object('down')
#ssort.navigate_selected_object('up')

###class keyev(object):
###    pass
###
###ev = keyev
###ev.key = 'down'
###ssort.key_on(ev)


#halt = True
