#%%
from multiprocessing.resource_sharer import stop
import sys
sys.path.append('sequence_jacobian')


import time


import importlib
import durtest
import household_models
importlib.reload(durtest)
importlib.reload(household_models)

from sequence_jacobian import create_model
from durtest import household_v, income, make_grids, firm, usercost, mkt_clearing, firm_ss
from household_models import household_nofc, household_localnofc


#   
r = 0.01
amax = 200
dmax = 12
xmax = amax + dmax
common_params = {'eis': 1, 'delta': 0.05, 'alpha': 0.11, 'rhoY': 0.966, 'sigmaY': 0.5, 'sigmaV': 0.00025,
                   'psi': 0.9, 'deltad': 1 - 0.8**(1/4), 'xi': 1, 'r': r, 'spread': 0.02, 'probadj': 0.99,
                   'f': 0, 'op_cost': 0.055, 'maint_cost': 0.466, 'collateral': 0.8, 'pbar': 1, 'pshock': 1,
                   'Y': 1.0, 'L': 1.0, 'nE': 2, 'nA': 100, 'amax': amax, 'nD': 50, 'dmax': dmax, 'dmin': 0.05, 'phiD': 1, 
                   'nX': 200, 'xmax': xmax, 'xmin': 0, 'nSfine': 40, 'nM': 200} 

def dagfc(vainit=0, vdinit=0, vinit=0):
    # Combine blocks
    household = household_v.add_hetinputs([income, make_grids])
    model = create_model([household, firm, usercost, mkt_clearing], name="FC Code")
    model_ss = create_model([household, firm_ss, usercost, mkt_clearing], name="FC Code")

    # Steady state
    
    calibration = {**common_params, 'vainit': vainit, 'vdinit': vdinit, 'vinit': vinit}
    unknowns_ss = {'beta': (0.985/(1+r), 0.9975/(1+r))}
    targets_ss = {'asset_mkt': 0.}
    ss = model_ss.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='brentq')


    return model_ss, ss, model

def dagnofc(vainit=0, vdinit=0, vinit=0):
    # Combine blocks
    household = household_nofc.add_hetinputs([income, make_grids])
    model = create_model([household, firm, usercost, mkt_clearing], name="FC Code")
    model_ss = create_model([household, firm_ss, usercost, mkt_clearing], name="FC Code")

    # Steady state
    calibration = {**common_params, 'vainit': vainit, 'vdinit': vdinit, 'vinit': vinit}
    unknowns_ss = {'beta': (0.985/(1+r), 0.9975/(1+r))}
    targets_ss = {'asset_mkt': 0.}
    ss = model_ss.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='brentq')


    return model_ss, ss, model  

def dagnofclocal(vainit=0, vdinit=0, vinit=0):
    # Combine blocks
    household = household_localnofc.add_hetinputs([income, make_grids])
    model = create_model([household, firm, usercost, mkt_clearing], name="FC Code")
    model_ss = create_model([household, firm_ss, usercost, mkt_clearing], name="FC Code")

    # Steady state
    calibration = {**common_params, 'vainit': vainit, 'vdinit': vdinit, 'vinit': vinit}
    unknowns_ss = {'beta': (0.975/(1+r), 0.9975/(1+r))}
    targets_ss = {'asset_mkt': 0.}
    ss = model_ss.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='brentq')

    return model_ss, ss, model        


if __name__ == '__main__':

    t0 = time.time()

    nofc_model_ss, nofc_ss, nofc_model = dagnofc()    

    # vainit = nofc_ss.internals['household_nofc']['Va']
    # vdinit = nofc_ss.internals['household_nofc']['Vd']
    # vinit = nofc_ss.internals['household_nofc']['Vd']

    nofclocal_model_ss, nofclocal_ss, nofclocal_model = dagnofclocal() #vainit, vdinit, vinit

    # vainit = nofclocal_ss.internals['household_localnofc']['Va']
    # vdinit = nofclocal_ss.internals['household_localnofc']['Vd']
    # vinit = nofclocal_ss.internals['household_localnofc']['Vd']

    fc_model_ss, fc_ss, fc_model = dagfc() #vainit, vdinit, vinit
    
    
    t1 = time.time()  
    print(t1-t0) 

    # c = {'nofc': nofc_ss.internals['household_nofc']['c'], 'nofclocal': nofclocal_ss.internals['household_localnofc']['c']}
    # d = {'nofc': nofc_ss.internals['household_nofc']['d1'], 'nofclocal': nofclocal_ss.internals['household_localnofc']['d1']} 

    # # assert abs(fc_ss['C'] - nofc_ss['C'])<0.01
    # # assert abs(fc_ss['D1'] - nofc_ss['D1'])<0.01

    # # assert abs(fc_ss['C'] - nofclocal_ss['C'])<0.01
    # # assert abs(fc_ss['D1'] - nofclocal_ss['D1'])<0.01

    # print(d['nofc'][0,0,0:5],d['nofclocal'][0,0,0:5])
    # print(d['nofc'][0,0,0:5]/c['nofc'][0,0,0:5],d['nofclocal'][0,0,0:5]/c['nofclocal'][0,0,0:5])
    # print(d['nofc'][-1,-1,0:5],d['nofclocal'][-1,-1,0:5])
    # print(d['nofc'][-1,-1,0:5]/c['nofc'][-1,-1,0:5],d['nofclocal'][-1,-1,0:5]/c['nofclocal'][-1,-1,0:5])

    # # Transitional dynamics
    # exogenous = ['Z', 'pshock']
    # unknowns = ['K']
    # targets = ['asset_mkt']

    # T = 200
    # nofc_G = nofc_model.solve_jacobian(nofc_ss, unknowns, targets, exogenous, T=T)
    # nofclocal_G = nofclocal_model.solve_jacobian(nofclocal_ss, unknowns, targets, exogenous, T=T)

    # print(nofc_G['D1']['pshock'][0:2,0:2],nofclocal_G['D1']['pshock'][0:2,0:2])


    # STOP


    c = {'fc': fc_ss.internals['household_v']['c'], 'nofc': nofc_ss.internals['household_nofc']['c'], 'nofclocal': nofclocal_ss.internals['household_localnofc']['c']}
    d = {'fc': fc_ss.internals['household_v']['d1'], 'nofc': nofc_ss.internals['household_nofc']['d1'], 'nofclocal': nofclocal_ss.internals['household_localnofc']['d1']} 

    # assert abs(fc_ss['C'] - nofc_ss['C'])<0.01
    # assert abs(fc_ss['D1'] - nofc_ss['D1'])<0.01

    # assert abs(fc_ss['C'] - nofclocal_ss['C'])<0.01
    # assert abs(fc_ss['D1'] - nofclocal_ss['D1'])<0.01

    print(d['fc'][0,0,0:5],d['nofc'][0,0,0:5],d['nofclocal'][0,0,0:5])
    print(d['fc'][-1,-1,0:5],d['nofc'][-1,-1,0:5],d['nofclocal'][-1,-1,0:5])
    

    # Transitional dynamics
    exogenous = ['Z', 'pshock']
    unknowns = ['K']
    targets = ['asset_mkt']

    T = 200
    fc_G = fc_model.solve_jacobian(fc_ss, unknowns, targets, exogenous, T=T)
    nofc_G = nofc_model.solve_jacobian(nofc_ss, unknowns, targets, exogenous, T=T)
    nofclocal_G = nofclocal_model.solve_jacobian(nofclocal_ss, unknowns, targets, exogenous, T=T)

    print(fc_G['D1']['pshock'][0:2,0:2],nofc_G['D1']['pshock'][0:2,0:2],nofclocal_G['D1']['pshock'][0:2,0:2])

    
# %%
