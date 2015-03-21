"""
This module holds the methods used to build a Model/Generator for a Hodgkin-
Huxley model of a neuron.  The user can build custom neurons by writing a config
file and importing it through the main run file.
"""

from __future__ import division
from numpy import *

from PyDSTool import *
from PyDSTool.Toolbox.dssrt import *
from PyDSTool.Toolbox.neuralcomp import *
from PyDSTool.Toolbox.phaseplane import find_fixedpoints, fixedpoint_2D, nullcline

from fovea.common import *

from PyDSTool import ModelManager

global man


def make_IN_model(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Type I excitability.
    """
    v = Var(voltage)
    ma = 0.32*(v+54.)/(1-Exp(-(v+54.)/4))
    mb = 0.28*(v+27.)/(Exp((v+27.)/5)-1)
    ha = .128*Exp(-(50.+v)/18)
    hb = 4/(1+Exp(-(v+27.)/5))
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', False, ma, mb, 3,
                                    'h', False, ha, hb, 1, vrev=50, g=100,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m','h'),
                                            'm': (voltage,),
                                            'h': (voltage,)})

    na = .032*(v+52)/(1-Exp(-(v+52.)/5))
    nb = .5*Exp(-(57.+v)/40)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4,
                                   vrev=-100,
                                   g=80, noauxs=False, subclass=channel_on,
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 2, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-67, g=0.1, noauxs=False,
                                  subclass=channel_on,
                                  gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'IN_typeI', targetGen,
                       with_spike_ev)


def make_IN_model_no_h(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Type I excitability with no h term.
    For spike onset only!
    """
    v = Var(voltage)
    ma = 0.32*(v+54.)/(1-Exp(-(v+54.)/4))
    mb = 0.28*(v+27.)/(Exp((v+27.)/5)-1)
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', False, ma, mb, 3,
                                    vrev=50, g=100,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m',),
                                            'm': (voltage,)})

    na = .032*(v+52)/(1-Exp(-(v+52.)/5))
    nb = .5*Exp(-(57.+v)/40)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4,
                                   vrev=-100,
                                   g=80, noauxs=False, subclass=channel_on,
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 2, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-67, g=0.1, noauxs=False,
                                  subclass=channel_on,
                                  gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'IN_typeI_no_h', targetGen,
                       with_spike_ev)


def make_HH_model(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Classic Type II excitability"""
    v = Var(voltage)
    ma = 0.1*(v+40)/(1-Exp(-(v+40)/10))
    mb = 4*Exp(-(v+65)/18)
    ha = .07*Exp(-(v+65)/20)
    hb = 1/(1+Exp(-(v+35)/10)) #0.01*(v+55)/(1-Exp(-(v+55)/10))
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', False, ma, mb, 3,
                                    'h', False, ha, hb, 1, vrev=50, g=120,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m','h'),
                                            'm': (voltage,),
                                            'h': (voltage,)})

    na = .01*(v+55)/(1-Exp(-(v+55)/10))
    nb = .125*Exp(-(v+65)/80)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-77,
                                   g=36, noauxs=False, subclass=channel_on,
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 8, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-54.4, g=0.3, noauxs=False,
                                    subclass=channel_on,
                                    gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'classic_typeII', targetGen,
                       with_spike_ev)


def make_HH_model_no_h(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Classic type II excitability with no h term in Na.
    Used for spike onsets only!"""
    v = Var(voltage)
    ma = 0.1*(v+40)/(1-Exp(-(v+40)/10))
    mb = 4*Exp(-(v+65)/18)
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', False, ma, mb, 3,
                                    vrev=50, g=120,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m',),
                                            'm': (voltage,)
                                            })

    na = .01*(v+55)/(1-Exp(-(v+55)/10))
    nb = .125*Exp(-(v+65)/80)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-77,
                                   g=36, noauxs=False, subclass=channel_on,
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 8, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-54.4, g=0.3, noauxs=False,
                                    subclass=channel_on,
                                    gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'classic_typeII_no_h', targetGen,
                       with_spike_ev)


def make_HH_model_fast_m_no_h(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Classic type II excitability with no h term in Na, fast m dynamics,
    and a constant time-scale K channel.
    Used for spike onsets only!"""
    v = Var(voltage)
    ma = 0.1*(v+40)/(1-Exp(-(v+40)/10))
    mb = 4*Exp(-(v+65)/18)
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', True, ma, mb, 3,
                                    vrev=50, g=120,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m',),
                                            'm': (voltage,)
                                            })

    na = .01*(v+55)/(1-Exp(-(v+55)/10))
    nb = .125*Exp(-(v+65)/80)
    ninf = na/(na+nb)
    channel_K1 = makeChannel_halfact('K', voltage, 'n', False, ninf, 'taun_par', 4, vrev=-77,
                                   g=36, noauxs=False, subclass=channel_on,
                                   parlist=[Par('5', 'taun_par')],
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 8, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-54.4, g=0.3, noauxs=False,
                                    subclass=channel_on,
                                    gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'classic_typeII_fast_m_no_h', targetGen,
                       with_spike_ev)




def make_HH_morepars_model_no_h(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Classic type II excitability with no h term in Na.
    Used to spike onsets only!"""
    v = Var(voltage)
    ma = 0.1*(v+40)/(1-Exp(-(v+40)/10))
    mb = 4*Exp(-(v+65)/18)
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', False, ma, mb, 3,
                                    vrev=50, g=120,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m',),
                                            'm': (voltage,)
                                            })

    ca = Par('0.01', 'ca')
    cb = Par('0.125', 'cb')
    aq = Par('10', 'aq')
    bq = Par('80', 'bq')
    va = Par('-55', 'va')
    vb = Par('-65', 'vb')
    na = ca*(v-va)/(1-Exp((va-v)/aq))
    nb = cb*Exp((vb-v)/bq)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-77,
                                   g=36, noauxs=False, subclass=channel_on,
                                   parlist=[ca,cb,aq,bq,va,vb],
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 8, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-54.4, g=0.3, noauxs=False,
                                    subclass=channel_on,
                                    gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'classic_typeII_morepars_no_h', targetGen,
                       with_spike_ev)


def make_WB_model(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Wang-Buzsaki Type I interneuron model"""
    phi = 5
    v = Var(voltage)
    ma = 0.1*(v+35)/(1-Exp(-(v+35)/10))
    mb = 4*Exp(-(v+60)/18)
    ha = phi*.07*Exp(-(v+58)/20)
    hb = phi*1/(1+Exp(-(v+28)/10))
    # isinstant=True for m
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', True, ma, mb, 3,
                                    'h', False, ha, hb, 1, vrev=55, g=35,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m','h'),
                                            'm': (voltage,),
                                            'h': (voltage,)})

    na = phi*.01*(v+34)/(1-Exp(-(v+34)/10))
    nb = phi*.125*Exp(-(v+44)/80)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-90,
                                   g=9, noauxs=False, subclass=channel_on,
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 2.5, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-65, g=0.1, noauxs=False,
                                    subclass=channel_on,
                                    gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'WB_typeI', targetGen,
                       with_spike_ev)


def make_WB_model_no_h(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Wang-Buzsaki Type I interneuron model with no h term in Na.
    Used to spike onsets only!"""
    phi = 5
    v = Var(voltage)
    ma = 0.1*(v+35)/(1-Exp(-(v+35)/10))
    mb = 4*Exp(-(v+60)/18)
    # isinstant=True for m
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', True, ma, mb, 3,
                                    vrev=55, g=35,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m',),
                                            'm': (voltage,)})

    na = phi*.01*(v+34)/(1-Exp(-(v+34)/10))
    nb = phi*.125*Exp(-(v+44)/80)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-90,
                                   g=9, noauxs=False, subclass=channel_on,
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 2.5, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-65, g=0.1, noauxs=False,
                                    subclass=channel_on,
                                    gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'WB_typeI_no_h', targetGen,
                       with_spike_ev)


def make_WB_model_minflin_no_h(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Wang-Buzsaki Type I interneuron model with no h term in Na,
    and a linear minf(V).
    Used for spike onsets only!"""
    phi = 5
    v = Var(voltage)
    v0 = -65.1
    m0 = 0.0285
    slope = 0.0035
    minf_quant = slope*(v-v0)+m0
    minf = 'max(0, %s)' % str(minf_quant)
    channel_Na1 = makeChannel_halfact('Na', voltage, 'm', True, minf, spow=3,
                                    vrev=55, g=35,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m',),
                                            'm': (voltage,)})

    na = phi*.01*(v+34)/(1-Exp(-(v+34)/10))
    nb = phi*.125*Exp(-(v+44)/80)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-90,
                                   g=9, noauxs=False, subclass=channel_on,
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 2.5, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-65, g=0.1, noauxs=False,
                                    subclass=channel_on,
                                    gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'WB_typeI_minflin_no_h', targetGen,
                       with_spike_ev)


def make_IN_morepars_model(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Type I HH model with kinetics pars for K channel"""
    v = Var(voltage)
    ca = Par('0.032', 'ca')
    cb = Par('0.5', 'cb')
    aq = Par('5', 'aq')
    bq = Par('40', 'bq')
    va = Par('-52', 'va')
    vb = Par('-57', 'vb')
    sw_h = Par('1', 'sw_h')
    ma = 0.32*(v+54.)/(1-Exp(-(v+54.)/4))
    mb = 0.28*(v+27.)/(Exp((v+27.)/5)-1)
    ha = sw_h*.128*Exp(-(50.+v)/18)
    hb = 4/(1+Exp(-(v+27.)/5))
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', False, ma, mb, 3,
                                    'h', False, ha, hb, 1, vrev=50, g=100,
                                    noauxs=False, subclass=channel_on,
                                    parlist=[sw_h],
                                    gamma1={voltage: ('m','h'),
                                            'm': (voltage,),
                                            'h': (voltage,)})

    na = ca*(v-va)/(1-Exp((va-v)/aq))
    nb = cb*Exp((v-vb)/bq)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-100,
                                   g=80, noauxs=False, subclass=channel_on,
                                   parlist=[ca,cb,aq,bq,va,vb],
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 2.1, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-67, g=0.1, noauxs=False,
                                    subclass=channel_on, gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'typeI_morepars', targetGen, with_spike_ev)


def make_IN_morepars_model_no_h(ic_args, alg_args, targetGen, with_spike_ev=False):
    """Type I HH model with kinetics pars for K channel.
    h removed to study spike onset.
    """
    v = Var(voltage)
    ca = Par('0.032', 'ca')
    cb = Par('0.5', 'cb')
    aq = Par('5', 'aq')
    bq = Par('40', 'bq')
    va = Par('-52', 'va')
    vb = Par('-57', 'vb')
    ma = 0.32*(v+54.)/(1-Exp(-(v+54.)/4))
    mb = 0.28*(v+27.)/(Exp((v+27.)/5)-1)
    channel_Na1 = makeChannel_rates('Na', voltage, 'm', False, ma, mb, 3,
                                    vrev=50, g=100,
                                    noauxs=False, subclass=channel_on,
                                    gamma1={voltage: ('m',),
                                            'm': (voltage,)})

    na = ca*(v-va)/(1-Exp((va-v)/aq))
    nb = cb*Exp((vb-v)/bq)
    channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-100,
                                   g=80, noauxs=False, subclass=channel_on,
                                   parlist=[ca,cb,aq,bq,va,vb],
                                   gamma1={voltage: ('n',),
                                           'n': (voltage,)})

    channel_Ib1 = makeBiasChannel('Ib', 2.1, noauxs=False, subclass=channel_on,
                                  gamma2={voltage: ('Ibias',)})
    channel_Lk1 = makeChannel_rates('Lk', vrev=-67, g=0.1, noauxs=False,
                                    subclass=channel_on, gamma1={voltage: ('Lk',)})

    return instantiate(ic_args, alg_args, [channel_Lk1, channel_Ib1,
                       channel_Na1, channel_K1], 'typeI_morepars_no_h',
                       targetGen, with_spike_ev)


make_model = {'HH_classic_typeII': make_HH_model,
 'HH_classic_typeII_no_h': make_HH_model_no_h,
 'HH_classic_typeII_fast_m_no_h': make_HH_model_fast_m_no_h,
 'HH_classic_typeII_morepars_no_h': make_HH_morepars_model_no_h,
 'HH_WB_typeI': make_WB_model,
 'HH_WB_typeI_no_h': make_WB_model_no_h,
 'HH_WB_typeI_minflin_no_h': make_WB_model_minflin_no_h,
 'HH_IN_typeI': make_IN_model,
 'HH_IN_typeI_no_h': make_IN_model_no_h,
 'HH_IN_typeI_morepars': make_IN_morepars_model,
 'HH_IN_typeI_morepars_no_h': make_IN_morepars_model_no_h}


def instantiate(ic_args, alg_args, channel_list, name, targetGen, with_spike_ev=False,
                withJac=False):
    """Presently, cannot use the Jacobian functionality with instantaneous activations.
    """
    global man
    man = ModelManager('HH_nullc_proj')

    if targetGen == 'Vode_ODEsystem':
        targetlang='python'
    else:
        targetlang='c'

    soma_name = 'cell_'+name

    HHcell = makeSoma(soma_name, channelList=channel_list, C=1.0, noauxs=False,
                            channelclass=channel_on)

    desc = GDescriptor(modelspec=HHcell, description='HH cell',
                       target=targetGen, algparams=alg_args)
    test = desc.validate()
    if not test[0]:
        print test[1]
        raise AssertionError
    assert desc.isinstantiable()

    # build an event that picks out when RHS of cell1's Voltage eqn is 0
    # i.e. when dV/dt=0, among others
    if with_spike_ev:
        HHcell.add(Par(0, 'V_th'))
        HHcell.flattenSpec()

        max_ev_args = {'name': 'cell_max',
                        'eventtol': 1e-8,
                        'eventdelay': 1e-3,
                        'starttime': 0,
                        'term': False
                        }
        # stationary event => dv/dt = 0
        max_ev = Events.makeZeroCrossEvent(HHcell.flatSpec['vars']['V'],
                                -1, max_ev_args, targetlang=targetlang,
                                flatspec=HHcell.flatSpec)
        min_ev = copy(max_ev)
        min_ev.dircode = 1
        min_ev.name = 'cell_min'

        # voltage threshold crossing event
        v_ev_args = {'name': 'V_thresh',
                        'eventtol': 1e-8,
                        'eventdelay': 1e-3,
                        'starttime': 0,
                        'term': False
                        }
        v_ev = Events.makeZeroCrossEvent('V-V_th', 0, v_ev_args,
                                         targetlang=targetlang,
                                         flatspec=HHcell.flatSpec)
        desc.userEvents = [min_ev, max_ev, v_ev]

    modelname = 'HH_'+name
    model_desc = MDescriptor(name=modelname,
                             withJac={HHcell.name: withJac})
    model_desc.add(desc)
    man.add(model_desc)

    # instantiate model for target generator
    # (could also try building syn_model_new)
    man.build(modelname, icvalues=ic_args, tdata=[0, 100])
    return man


def get_model(path, name, ic_args, alg_args, gen='Vode'):
    global man
    try:
        man = loadObjects(path+'model_'+name+'.sav')[0]
    except:
        # make_model is a dict for looking up modelspec factory functions,
        # ultimately calling instantiate and returning model manager 'man'
        man = make_model[name](ic_args, alg_args, gen+'_ODEsystem', with_spike_ev=False)
        saveObjects(man, path+'model_'+name+'.sav')
    return man


def get_ref(HHmodel, path, name, force=False):
    try:
        if force:
            raise "recalculate"
        ref_ic, ref_pts, ref_traj = loadObjects(path+'ref_'+name+'.sav')
    except:
        HHmodel.compute(trajname='ref_traj', tdata=[0,400])
        evts = HHmodel.getTrajEventTimes('ref_traj', 'cell_min')
        d_evts = [evts[i]-evts[i-1] for i in range(1,len(evts))]
        ref_ic = HHmodel('ref_traj', evts[-1])
        ref_pts = HHmodel.sample('ref_traj', tlo=evts[-2], thi=evts[-1])
        ref_pts.indepvararray -= evts[-2]
        ref_traj = pointset_to_traj(ref_pts)
        saveObjects([ref_ic, ref_pts, ref_traj], path+'ref_'+name+'.sav', force=force)
    return ref_ic, ref_pts, ref_traj



def getHH_DSSRT(model, pts, name=None, header=None, verbose=0,
                sigma=3, gamma=3, force=False):
    """
    DSSRT quantities (psi's and omega's) are defined by hand.
    After calculating quantities the points are returned as a pointset.

    Set force option True to avoid trying to load anything (e.g. in case ICs changed)
    """

    if header is None:
        # Select choice for contents of "op" post-processed pointset and order & choice of columns in string table
        header = ['t', 'Psi_l', 'Psi_i', 'Psi_m', 'Psi_n', 'Omega_m', 'Omega_n', 'Vdot', 'mdot', 'ndot', 'd_real',
          'A', 'accn', 'd', 'a', 'theta', 'c',
          'tleft', 'tleft_real', 'horiz_t', 'spike', 'tauv', 'taum', 'taun', 'V', 'm', 'n']


    if name is None:
        assert model.name in make_model.keys(), "Generator must be of pre-defined HH-model."
        name = model.name

    if model.name[-5:] == '_no_h':
        no_h = True
    else:
        no_h = False


    if verbose >= 1:
        print "Calculating DSSRT for: %s" % name

    gen = model.registry.values()[0]

    ############################

    if force:
        loaded_ok = False
    else:
        try:
            dssrt_data_tuple = loadObjects(name+'_DSSRT.sav')
        except:
            loaded_ok = False
        else:
            try:
                loaded_ok = all(dssrt_data_tuple[1].coordnames == \
                    remain(list(sort(header)), ['t'])) and \
                    all(dssrt_data_tuple[0] == pts)
            except:
                print "Re-run with force option True"
                raise


    # Can't save and reload dssrt_assistant objects yet, so re-create those

    # set up quantities for analysis (mainly DSSRT stuff)

    mspec = model._mspecdict[gen.name].modelspec

    Dargs = args(model=gen, use_order=True, speed_fussy=True)
    Dargs.inputs = mspec.dssrt_inputs
    # inputs (for V) corresponds to:
    #{'V': args(gamma1=['Na.m', 'Na.h', 'K.n', 'Lk.Lk'],  # unless Na.h removed
    #                     gamma2=['Ib.Ibias'])}
    Dargs.taus = {'V': 'tauv', 'Na.m': 'Na.taum', 'Na.h': 'Na.tauh',
                  'K.n': 'K.taun', 'Lk.Lk': None}
    Dargs.infs = {'V': 'vinf', 'Na.m': 'Na.minf', 'Na.h': 'Na.hinf',
                  'K.n': 'K.ninf', 'Lk.Lk': None}

    if no_h:
        del Dargs.infs['Na.h']
        del Dargs.taus['Na.h']

    # set up terms for dm/dt and Psi_m depending on whether m = m_inf(V) in the Wang-Buzsaki model
    if model.name[:11] == "HH_WB_typeI" or 'fast_m' in model.name:
        if 'minflin' in model.name:
            # linearized testing of m_inf
            mdot_string = '(1/tauv)*(vinf-V)*0.0035'
        else:
            mdot_string = '(1/tauv)*(vinf-V)*((4.0 + 0.1*V)*(-0.1/(1 - exp(-4 - V/10)) + (4.0 + 0.1*V)*exp(-4 - V/10)/(10*(1 - exp(-4 - V/10))**2) + 2*exp(-65./18 - V/18)/9)/((1 - exp(-4 - V/10))*((4.0 + 0.1*V)/(1 - exp(-4 - V/10)) + 4*exp(-65./18 - V/18))**2) + 0.1/((1 - exp(-4 - V/10))*((4.0 + 0.1*V)/(1 - exp(-4 - V/10)) + 4*exp(-65./18 - V/18))) - (4.0 + 0.1*V)*exp(-4 - V/10)/(10*(1 - exp(-4 - V/10))**2*((4.0 + 0.1*V)/(1 - exp(-4 - V/10)) + 4*exp(-65./18 - V/18))))'
        # the tauv's cancel here
        mterm = 'tauv*Na.g*3*Na.m*Na.m*Na.h*abs((Na.vrev-vinf)*' + mdot_string + ')'
    else:
        mterm = 'tauv*Na.g*3*Na.m*Na.m*Na.h*abs((Na.vrev-vinf)*(Na.minf-Na.m)/Na.taum)'
    if no_h:
        mterm = mterm.replace('*Na.h','')


    orig_omegas = {'V':
                {'Na.m': mterm,
                 'Na.h': 'tauv*Na.g*Na.m*Na.m*Na.m*abs((Na.vrev-vinf)*(Na.hinf-Na.h)/Na.tauh)',
                 'K.n': 'tauv*K.g*4*K.n*K.n*K.n*abs((K.vrev-vinf)*(K.ninf-K.n)/K.taun)'
                 }}

    orig_psis = {'V':
                {'Na.m': 'tauv*Na.g*3*Na.m*Na.m*Na.h*abs(Na.vrev-vinf)',
                 'Na.h': 'tauv*Na.g*Na.m*Na.m*Na.m*abs(Na.vrev-vinf)',
                 'K.n': 'tauv*K.g*4*K.n*K.n*K.n*abs(K.vrev-vinf)',
                 'Lk.Lk': 'tauv*Lk.g*abs(Lk.vrev-vinf)',
                 'Ib.Ibias': 'tauv*abs(Ib.Ibias)'
                }}

    if no_h:
        del orig_omegas['V']['Na.h']
        del orig_psis['V']['Na.h']
        orig_omegas['V']['Na.m'] = orig_omegas['V']['Na.m'].replace('*Na.h','')
        orig_psis['V']['Na.m'] = orig_psis['V']['Na.m'].replace('*Na.h','')


    # creates DSSRT objects that know how to compute Psi and Omega values
    traj = pointset_to_traj(pts)

    # unsigned influences (don't need here)
    Dargs.psis = orig_psis
    da_reg = dssrt_assistant(Dargs)
    da_reg.focus_var = 'V'
    da_reg.traj = traj

    Dargs.psis = orig_omegas
    da_rate = dssrt_assistant(Dargs)
    da_rate.focus_var = 'V'
    da_rate.traj = traj

    # signed influences
    Dargs_reg_signed = copy(Dargs)
    Dargs_reg_signed.psis = copy(orig_psis)
    for xinput, psi_def in Dargs_reg_signed.psis['V'].items():
        Dargs_reg_signed.psis['V'][xinput] = psi_def.replace('abs', '')
    da_reg_signed = dssrt_assistant(Dargs_reg_signed)
    da_reg_signed.focus_var = 'V'

    Dargs_rate_signed = copy(Dargs)
    Dargs_rate_signed.psis = copy(orig_omegas)
    for xinput, psi_def in Dargs_rate_signed.psis['V'].items():
        Dargs_rate_signed.psis['V'][xinput] = psi_def.replace('abs', '')
    da_rate_signed = dssrt_assistant(Dargs_rate_signed)
    da_rate_signed.focus_var = 'V'

    # ----------------------------------
    if loaded_ok:
        # (pts, ptsDSSRT, strDSSRT, epochs_reg, epochs_rate), (da_reg, da_rate, da_reg_signed, da_rate_signed)
        return dssrt_data_tuple[1:], (da_reg, da_rate, da_reg_signed, da_rate_signed)

    # else re-build everything
    da_reg.calc_psis()
    da_reg.make_pointsets()
    da_reg.calc_rankings()
    da_reg.domscales['psi'].calc_epochs(sigma, gamma)
    epochs_reg = da_reg.domscales['psi'].epochs

    da_rate.calc_psis()
    da_rate.make_pointsets()
    da_rate.calc_rankings()
    da_rate.domscales['psi'].calc_epochs(sigma, gamma)
    epochs_rate = da_rate.domscales['psi'].epochs

    # -----------------------

    # runtime-evaluated definitions for use in make_generic_table()
    # pairs are (defn_str, number_sig_figs) for string formatted output
    table_lookup = dict(
        t=("pts.indepvararray[i]", 3),
        Omega_m=("omega_m", 3),
        Omega_n=("omega_n", 3),
        Psi_m=("psi_m", 3),
        Psi_n=("psi_n", 3),
        Psi_l=("psi_l", 3),
        Psi_i=("psi_i", 3),
        Vdot=("rhs['V']", 4),
        ndot=("rhs['K.n']", 4),
        # mdot is later

        # change A when re-introduce Na.h (+ omega_h)
        A=("omega_m - rhs['V'] + omega_n", 3),

        tauv=("tauv", 3),
        taum=("taum", 3),
        taun=("taun", 3),
        vinf=("vinf", 3),
        minf=("minf", 3),
        ninf=("ninf", 3),

        # angles in linearized nullcline geometry derivation
        theta=("theta", 4),
        a=("a", 4),   # angle a already defined in make_generic_table
        c=("theta - a", 4),

        # distance to nullcline collision assuming linearized velocities and linearized nullcline
        d=("(vinf-pt['V'])*(sin(theta)/sin(a))", 3),
        dd_dt=("dd_dt", 4),

        # time to collision
        tleft=("d/dd_dt", 3),

        # projected "real" distance to nullcline collision using linearized velocities only
        d_real=("d_real", 3),
        tleft_real=("d_real/dd_dt", 3),

        # simply forward differencing of Vdot to find V's acceleration in time
        accn=("fdiff('Vdot')", 4),

        V=("pt['V']", 3),
        m=("pt['Na.m']", 3),
        n=("pt['K.n']", 3),

        # predictor for when A will become positive (only helpful when it is currently < 0!)
        horiz_t=("-_vals_['A']/fdiff('A')", 3),

        # spike (Boolean) is a predictor of whether a spike is imminent or not
        # it doesn't know if nullclines have crossed (when it fails to be valid)!
        spike=("int( (_vals_['horiz_t'] > 0 and _vals_['tleft'] > 0 and " + \
                       "_vals_['horiz_t'] < _vals_['tleft'] and fdiff('A') > 0 )    or " + \
                       "( fdiff('A') > 0 and _vals_['A']>0 and abs(psi_m/psi_n) > 2 and rhs['V'] > 0))", 0)
        )

    if name[:11] == "HH_WB_typeI" or 'fast_m' in name:
        table_lookup['mdot'] = (mdot_string.replace('V',"pt['V']"), 4)
    else:
        table_lookup['mdot'] = ("rhs['Na.m']", 4)

    ptsDSSRT, strDSSRT = make_generic_table(header, pts, model, da_reg_signed, da_rate_signed, table_lookup, repeat_header=65, verbose=verbose-1)
    saveObjects([pts, ptsDSSRT, strDSSRT, epochs_reg, epochs_rate], name+'_DSSRT.sav', force=True)
    print "\n"
    return (ptsDSSRT, strDSSRT, epochs_reg, epochs_rate), (da_reg, da_rate, da_reg_signed, da_rate_signed)


def make_generic_table(header, pts, model, da_reg, da_rate, table_lookup_dict, repeat_header=None, verbose=0):
    """If repeat_header is an integer, the header will repeat every that many rows
    """

    lenPts = len(pts)
    if verbose >= 1:
        print "Generating DSSRT points and table: %i pts" % lenPts

    str_table = []
    str_table.append(header)
    # hack to avoid using if statements in loop
    if repeat_header is None:
        repeat_header = NaN
    output_pts = {}
    for h in header:
        output_pts[h] = []

    # utility function for finite-differencing
    if 't' in header:
        def fdiff(h):
            vals = output_pts[h]
            try:
                dv = vals[-1] - vals[-2]
            except IndexError:
                # first point!
                return np.NaN
            else:
                ts = output_pts['t']
                return dv/(ts[-1]-ts[-2])
    else:
        def fdiff(h):
            vals = output_pts[h]
            try:
                dv = vals[-1] - vals[-2]
            except:
                # first point!
                return np.NaN
            else:
                return dv

    for i, pt in enumerate(pts):
        if verbose >= 1:
            progressBar(i, lenPts)
        if i > 0 and mod(i, repeat_header) == 0:
            str_table.append(header)
        row = []
        rhs = model.Rhs(0, pt, model.pars)
        psi_m = da_reg.calc_psi('Na.m', pt)
        psi_n = da_reg.calc_psi('K.n', pt)
        psi_l = da_reg.calc_psi('Lk.Lk', pt)
        psi_i = da_reg.calc_psi('Ib.Ibias', pt)
        omega_m = da_rate.calc_psi('Na.m', pt)
        omega_n = da_rate.calc_psi('K.n', pt)
        tauv = da_rate.calc_tau('V', pt)
        taum = da_rate.calc_tau('Na.m', pt)
        taun = da_rate.calc_tau('K.n', pt)
        vinf = da_rate.calc_inf('V', pt)
        minf = da_rate.calc_inf('Na.m', pt)
        ninf = da_rate.calc_inf('K.n', pt)
        theta = arctan2(1, psi_m)
        if 'Na.m' in rhs:
            a = theta - arctan2( rhs['Na.m'], (rhs['V']-omega_n) )
            dd_dt = sqrt((rhs['V']-omega_n)**2 + (rhs['Na.m'])**2)
        else:
            mdot = table_lookup_dict['mdot']
            a = theta - arctan2( eval(mdot[0]), (rhs['V']-omega_n) )
            dd_dt = sqrt((rhs['V']-omega_n)**2 + eval(mdot[0])**2)
        d = (vinf-pt['V'])*(sin(theta)/sin(a))

        mf = make_mline(pt, theta-a)
        V = pt['V']

        if 'd_real' in header:
            v_step = 0.25
            vhi = V
            vhi_max = 60
            if 'Na.h' in model.query('vars'):
                def obj_fn(V):
                    return da_reg.calc_inf('V', {'V': V, 'Na.m': mf(V), 'K.n': pt['K.n'], 'Na.h': pt['Na.h']}) - V
            else:
                def obj_fn(V):
                    return da_reg.calc_inf('V', {'V': V, 'Na.m': mf(V), 'K.n': pt['K.n']}) - V
            obj_a = vinf - V
            V_contact = None
            while vhi < vhi_max:
                vhi = vhi + v_step
                if obj_fn(vhi)*obj_a < 0:
                    V_contact = bisection(obj_fn, V, vhi, xtol=1e-4)
                    break

            if V_contact is None:
                d_real = np.NaN
            else:
                d_real = sqrt((V_contact-V)**2 + (mf(V_contact)-pt['Na.m'])**2)

        # convenience for reference by user
        _vals_ = {}
        for h in header:
            ev_str, sigd = table_lookup_dict[h]
            val = eval(ev_str)
            _vals_[h] = val
            output_pts[h].append(val)
            try:
                row.append(n_sigdigs_str(val, sigd))
            except:
                print "Error processing string %s for header %s" % (ev_str, h)
                raise

        str_table.append(row)
    tabulated = indent(str_table, hasHeader=True)
    if 't' in header:
        return Pointset(coorddict=filteredDict(output_pts, ['t'], neg=True),
                        indepvarname='t',
                        indepvararray=output_pts['t']), tabulated
    else:
        return Pointset(coorddict=output_pts), tabulated


from time import clock

def computePPlaneObjects(gen, xvar, yvar, state=None, tol=1e-8,
              x_scale=1, max_step=1, subdomain=None, seed_pts=None,
              only_var=None, do_fps=True, num_x_points=100, num_y_points=100,
              strict_domains=True, with_jac=True, fast_vars=None, verbose=0,
              as_array=True):
    """
    THIS VERSION IS HARDWIRED FOR USE WITH HODGKIN-HUXLEY MODELS.

    Compute basic phase plane objects (fixed points and nullclines)
    given certain parameters (system params will come from generator):

    gen         Generator object from which to build nullcline from
    xvar        Variable to be treated as independent
    yvar        Variable to be treated as dependent
    fast_vars   Any fast variables set to their adiabatic equilibrium values
    state       Point or dictionary of state values
    tol         Tolerance scalar
    max_step    Maximum step in arc-length that nullcline computation
                will use
    subdomain   Optional sub-domain of calculation (default to whole
                variable domains given by Generator) given as dictionary
                with [xlo, xhi] values for each x string key.
    seed_pts    Seed points (list of dictionaries or Points)
                for any known nullcline points on either/both nullclines.
    only_var    Select one variable to compute nullclines for only.
    do_fps      Boolean switch (default True) to control computation of
                fixed points.
    num_x_points    number of mesh points for x nullcline (default 100)
    num_y_points    number of mesh points for y nullcline (default 100)
    strict_domains  Boolean switch (default True) to control strict domain
                    cropping of nullclines computed.
    with_jac    Toggle to create Jacobian for use in fixed point and nullcline
                calculations (default True)
    verbose     0 = No text output; No graphic output
                1 = Text output only
                2 = Graphic output from this function only
                    No grahphic output from functions called within
                3+ = Pull output from subsequently lower-level functions
                    called within
    as_array    Optional Boolean to return nullclines as array data rather
                than fully interpolated spline nullcline objects (default True)
    """

    ## ------------------------------------------------------------
    ## JACOBIAN
    ## ------------------------------------------------------------

    state = dict(state)  # ensure updatable!
    yvarFS = gen._FScompatibleNames(yvar)
    xvarFS = gen._FScompatibleNames(xvar)
    stateFS = filteredDict(gen._FScompatibleNames(state), gen.funcspec.vars)
    varPars = filteredDict(stateFS, [xvarFS, yvarFS], neg=True)

    # TEMP
    #with_jac = False
    #do_fps = False

    if with_jac:
        jacFS, new_fnspecsFS = prepJacobian(gen.funcspec._initargs['varspecs'], [xvarFS, yvarFS],
                                    gen.funcspec._initargs['fnspecs'])

        scopeFS = gen.pars.copy()
        scopeFS.update(varPars)
        scopeFS.update(new_fnspecsFS)
        jac_fnFS = expr2fun(jacFS, ensure_args=['t'], **scopeFS)

        jac = gen._FScompatibleNamesInv(jacFS)
        scope = gen._FScompatibleNamesInv(gen.pars)
        scope.update(filteredDict(state, [xvar, yvar], neg=True))
        new_fnspecs = {}
        # don't convert the new_fnspecs keys to . notation
        keymap = Symbolic.symbolMapClass(dict(zip(gen._FScompatibleNamesInv(new_fnspecsFS.keys()), new_fnspecsFS.keys())))
        for k, v in new_fnspecsFS.items():
            new_fnspecs[k] = Symbolic.mapNames(keymap, gen._FScompatibleNamesInv(v))
        scope.update(new_fnspecs)
        jac_fn = expr2fun(jac, ensure_args=['t'], **scope)
    else:
        jac_fn = jac_fnFS = None

    ## ------------------------------------------------------------
    ## FIXED POINTS
    ## ------------------------------------------------------------

    if verbose >= 1:
        print "Computing fixed points"

    if subdomain is not None:
        stateFS.update(filteredDict(gen._FScompatibleNames(subdomain),
                                    (xvarFS, yvarFS)))

    if do_fps:
        # create Intervals to test inclusion in subdomain
        int_x = Interval(xvar, float, subdomain[xvar])
        int_y = Interval(yvar, float, subdomain[yvar])

        # temporarily overwriting state with xdomains from gen,
        # for fussy fsolve in find_fixedpoints
        for v in (xvarFS, yvarFS):
            stateFS[v] = gen.xdomain[v]
        fp_coords = find_fixedpoints(gen, n=10, jac=jac_fnFS, subdomain=stateFS,
                                 eps=tol)

        fps = []
        for fp in fp_coords:
            if fp[xvar] in int_x and fp[yvar] in int_y:
                fps.append(fixedpoint_2D(gen, Point(fp), coords=[xvar, yvar],
                                 jac=jac_fn, description='', eps=1e-5))
        # reset stateFS to use subdomain
        if subdomain is not None:
            stateFS.update(filteredDict(gen._FScompatibleNames(subdomain),
                                    (xvarFS, yvarFS)))
    else:
        fp_coords = None
        fps = []


    ## ------------------------------------------------------------
    ## NULLCLINES
    ## ------------------------------------------------------------

    gen.set(ics=varPars)

    if verbose >= 1:
        print "Computing nullclines"

    assert 'V' in (xvar, yvar)
    if 'V' == xvar:
        #Vs = linspace(gen.xdomain['V'][0], gen.xdomain['V'][1], num_x_points)
        Vs = linspace(stateFS['V'][0], stateFS['V'][1], num_x_points)
        other = yvar
        otherFS = yvarFS
        os = linspace(stateFS[yvarFS][0], stateFS[yvarFS][1], num_y_points)
        n = num_y_points
    else:
        Vs = linspace(stateFS['V'][0], stateFS['V'][1], num_y_points)
        other = xvar
        otherFS = xvarFS
        os = linspace(stateFS[xvarFS][0], stateFS[xvarFS][1], num_x_points)
        n = num_x_points

    # yes, this gets recalc'd even if only_var is 'V' but this is a cheap calc compared to
    # nullcline object creation
    ofn = getattr(gen.auxfns, other.split('.')[0]+'_dssrt_fn_'+other.split('.')[1]+'inf')
    oinfs = [ofn(V) for V in Vs]

    fn_args = gen.funcspec._initargs['fnspecs']['dssrt_fn_Vinf'][0]
    # Dictionary to map Na_m into m, etc.
    model_vars = gen.query('vars')
    varnamemap = Symbolic.symbolMapClass(dict([(v,v.split('.')[-1]) \
                                               for v in model_vars]))
    pt = filteredDict(Symbolic.mapNames(varnamemap,
                                        filteredDict(state, model_vars)),
                      fn_args)
    vinfs = np.zeros((n,),float)

    oarg = varnamemap[other]
    if other not in model_vars and fast_vars is not None and other in fast_vars:
        oarg = oarg.split('.')[-1]
        os = oinfs

    if fast_vars is not None and 'Na.m' in fast_vars:
        pt['m'] = gen.auxfns.Na_dssrt_fn_minf(state['V'])

    for i, o in enumerate(os):
        pt[oarg] = o
        vinfs[i] = gen.auxfns.dssrt_fn_Vinf(**pt)

    if 'V' == xvar:
        nulls_x = array([vinfs, os]).T
        nulls_y = array([Vs, oinfs]).T
    else:
        nulls_y = array([os, vinfs]).T
        nulls_x = array([oinfs, Vs]).T

    if as_array:
        nullcX = nulls_x
        nullcY = nulls_y
    else:
        if only_var is None:
            nullcX = nullcline(xvar, yvar, nulls_x, x_relative_scale=x_scale)
            nullcY = nullcline(xvar, yvar, nulls_y, x_relative_scale=x_scale)
        elif only_var == xvar:
            nullcX = nullcline(xvar, yvar, nulls_x, x_relative_scale=x_scale)
            nullcY = None
        elif only_var == yvar:
            nullcX = None
            nullcY = nullcline(xvar, yvar, nulls_y, x_relative_scale=x_scale)
        else:
            raise ValueError("Invalid variable name for only_var: %s" % only_var)

    return {'nullcX': nullcX, 'nullcY': nullcY, 'fps': fps}


def do_traj(model, t_end, trajname='test', cutoff_after_max=1, do_plot=True):
    """Utility function to compute a trajectory in the given model,
    optionally truncating the pointset (not the trajectory) after
    cutoff_after_max time units after a maximum occurs in voltage.

    Returns traj, pts (where pts is truncated to a range starting from
    increasing voltage to up to cutoff_after_max ms after a max)

    This function also plots the (maybe truncated) pointset in Figure 100.
    """
    model.compute(trajname, tdata=[0,t_end], force=True)
    traj = model[trajname]

    try:
        tmax = traj.eventTimes['cell_max'][0]
    except:
        pts = traj.sample()
        if t_end > 15 and pts[-1]['V'] < 0:
            print "Warning: expect monotonic approach to a sub-threshold f.p. in V, at least up to t=", t_end
            tmax = t_end
        else:
            raise RuntimeError("Did not run for long enough - V never reached a max")
    else:
        print "V had a max at t=",tmax, "with value", traj.events['cell_max'][0]['V']
    try:
        tmin = traj.eventTimes['cell_min'][0]
    except:
        tlo = 0
    else:
        if tmin < tmax:
            # ignore leading time where V is decreasing
            tlo = tmin
        else:
            tlo = 0
    if cutoff_after_max is not None:
        assert cutoff_after_max > 0
        pts = traj.sample(tlo=tlo, thi=tmax+cutoff_after_max)
    else:
        pts = traj.sample()
    if do_plot:
        pp.figure(100)
        pp.plot(pts['t'], pts['V'])
    return traj, pts, tlo, tmax
