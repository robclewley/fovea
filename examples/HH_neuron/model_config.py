"""
This is a template config file that can be used to build a Hodgkin-Huxley neuron
model.
"""

global man
from HH_neuron import get_model


## ----- ----- ----- ----- ----- ----- ##
## SET INITIAL CONDITIONS              ##
## ----- ----- ----- ----- ----- ----- ##

path = 'models/'


## ----- ----- ----- ----- ----- ----- ##
## SELECT ONE MODEL HERE               ##
## ----- ----- ----- ----- ----- ----- ##

name = 'HH_classic_typeII'
#name = 'HH_classic_typeII_no_h'
#name = 'HH_classic_typeII_fast_m_no_h'

#name = 'HH_classic_typeII_morepars_no_h'
#name = 'HH_IN_typeI_morepars_no_h' # no model for this yet
#name = 'HH_IN_typeI_no_h'

## Don't bother with WB for now!
#name = 'HH_WB_typeI_no_h'

#name = 'HH_IN_typeI'
#name = 'HH_classic_typeII'
#name = 'HH_WB_typeI'

## ----- ----- ----- ----- ----- ----- ##
## SET INITIAL CONDITIONS              ##
## ----- ----- ----- ----- ----- ----- ##

if 'typeI' in name and 'typeII' not in name:
    ic_args = {'V':-70.0, 'Na.m': 0,
           'Na.h': 0.7, 'K.n': 0.2}
else:
    ic_args = {'V':-70.0, 'Na.m': 0,
           'Na.h': 0.45, 'K.n': 0.46}


## ----- ----- ----- ----- ----- ----- ##
## SET SIMULATION PARAMETERS           ##
## ----- ----- ----- ----- ----- ----- ##

alg_args = {'init_step': 0.01, 'max_step': 0.05,
            'max_pts': 25000}


print("Using neuron model ", name)

# have to eliminate h from ICs if model does not use h
ic_args_no_h = ic_args.copy()
del ic_args_no_h['Na.h']

if name[-5:] == '_no_h':
    no_h = True
    man = get_model(path, name, ic_args_no_h, alg_args, gen='Dopri')
else:
    no_h = False
    man = get_model(path, name, ic_args, alg_args, gen='Dopri')

HHmodel = man.instances[name]
# ensure these are set (in case model was reloaded not created)
HHmodel.set(algparams=alg_args)
if no_h:
    HHmodel.set(ics=ic_args_no_h)
else:
    HHmodel.set(ics=ic_args)

if no_h:
    # h set to some constant => gNa adjusted
    h_const = ic_args['Na.h']
    orig_gNa = HHmodel.pars['Na.g']
    HHmodel.set(pars={'Na.g': orig_gNa*h_const})
    print("Na.g adjusted by h_const factor, original value stored in orig_gNa")

##### Header for selecting derived quantities to compute
# order can matter in this list, as later quantities may depend on earlier ones being up to date
# e.g. t needs to be before anything taking a finite difference
header = ['t', 'theta', 'a', 'c', 'd', 'dd_dt', 'A', 'Vdot', 'mdot', 'ndot', 'accn', 'tleft', 'horiz_t',
          'Psi_l', 'Psi_i', 'Psi_m', 'Psi_n', 'Omega_m', 'Omega_n']