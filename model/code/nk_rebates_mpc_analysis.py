#%%
import numpy as np
import pandas as pd
import pickle

# these are the sequence space packges
import utils
from simple_block import simple
import jacobian as jac
import nonlinear

# these are the model files
from nk_rebates_mpc_monthly_ss import nkrebate
from nk_rebates_mpc_model import investment, household_o, household_r, pricesetting, union, monetary, fiscal, marketclearing, durablesupply

# ------------------------------------------------------------------------
# Spending profile of tax shock
# ------------------------------------------------------------------------

# rebate scaled as a fraction of PCE
taxshock = pd.read_stata('../input/rel_rebate.dta')['rel_rebate'].to_numpy()

Tmax = 12
# taxshock = np.zeros(Tmax)
# taxshock[3] = 0.002
# taxshock[4] = 0.039
# taxshock[5] = 0.023
# taxshock[6] = 0.011
# taxshock[7] = 0.0008
# taxshock[8] = 0.0005
# taxshock[9] = 0.0008
# taxshock[10] = 0.0009

ratelist = ['r','pi','rr','rk']

# ------------------------------------------------------------------------
# Variable Parameters
# ------------------------------------------------------------------------

intertempelasticity = 15
intertempelasticitybkmm = 6.4
lrhorizon = 36
mpctargetsbase = [0.3, 0.5, 0.9]
mpcdurtargetsbase = [0.3, 0.4, 0.4]
mpcnonduradd = [0.06]

parasetnames = ['baseline','durableprice','lesselastic','nondurablesonly','sameexponent','sameexponentdurableprice','sameexponentdurablepricelow','sameexponentdurablepricehigh','inelasticdurables','roffset'] # #'gefeedback',,'inelasticsupply','strongmp' #['inelasticdurables']#
parasets = [{'zeta': 0},{'zeta': 0.2},{'zeta': 0.2},{'xc': 0},{'sigma': 1, 'sigmad': 1,'varpi': 1, 'xi':1}, {'sigma': 1, 'sigmad': 1,'varpi': 1, 'xi':1, 'zeta': 0.2}, {'sigma': 1, 'sigmad': 1,'varpi': 1, 'xi': 0.5, 'zeta': 0.2}, {'sigma': 1, 'sigmad': 1,'varpi': 1, 'xi': 1.5, 'zeta': 0.2},{'zeta': 1},{'phigap': 100, 'rhor': 0}] 


# parasetnames = ['inelasticdurables'] # #'gefeedback',,'inelasticsupply','strongmp' #['inelasticdurables']#
# parasets = [{'zeta': 1}] 

# parasetnames = ['roffset'] 
# parasets = [{'phigap': 100, 'rhor': 0}] 

# parasetnames = ['baseline'] 
# parasets = [{'zeta': 0}] 

parameter_set = ['gamma'] 
info_set = ['no-anticipation'] #['full-info'] 

model_names = dict()
model_names['Real'] = ['exp', 'realx', 'realc']
model_names['Nominal'] = ['nomexp', 'nomx', 'nomc']

# ------------------------------------------------------------------------
# Set-ups for solving the model
# ------------------------------------------------------------------------

T = {'micro': 60, 'GE': 200}

block_list = dict()
block_list['micro'] = [household_o, household_r, marketclearing, fiscal] #
block_list['GE'] = [investment, household_o, union, pricesetting, monetary, fiscal, marketclearing, household_r, durablesupply]

def getlists(block_list):
    allexogenous = set()
    allunknowns = set()
    targets = []
    for block in block_list:
        allexogenous = allexogenous | block.exogenous
        allunknowns = allunknowns | block.unknowns
        targets = targets + block.output_list

    exog = allexogenous - allunknowns
    unknowns = allunknowns - exog 
    unknowns = list( unknowns )
    exog = list( exog )

    return exog, unknowns, targets

exog=dict()
unknowns=dict()
targets=dict()

for key in T.keys():
    exog[key], unknowns[key], targets[key] = getlists(block_list[key])




# ------------------------------------------------------------------------
# Steady state parameters that match MPCs
# ------------------------------------------------------------------------


for parasetname, paraset in zip(parasetnames, parasets):
    print(paraset, parasetname)

    if parasetname=='nondurablesonly':
        # combine lists of mpc targets

        mpctargets = mpcnonduradd + mpctargetsbase
        mpcdurtargets = mpcnonduradd + mpcdurtargetsbase
    else:
        mpctargets = mpctargetsbase
        mpcdurtargets = mpcdurtargetsbase

    # find MPC
    mpc = dict()
    G = dict()
    output = dict()
    
    for parameter in parameter_set:

        thetad = 0.8
        
        for mpctarget, mpcdurtarget in zip(mpctargets, mpcdurtargets):
            print(parameter, mpctarget)

            gamma = mpctarget
            mpx = mpcdurtarget / mpctarget
            vartheta = 0
            fc = 0
            diff = 1

            maxit = 400

            for kk in range(maxit):

                ss = nkrebate(gamma=gamma, thetad=thetad, mpx=mpx, **paraset) 
                
                G['micro'] = jac.get_G(block_list=block_list['micro'],
                                exogenous=exog['micro'],
                                unknowns=unknowns['micro'],
                                targets=targets['micro'],
                                T=T['micro'], ss=ss)

                diffmpc  = G['micro']['exp']['incshock'][0:3,0].sum() - mpctarget
                if parasetname=='lesselastic':
                    diffinter = G['micro']['realx']['pdur'][0:2,2:T['micro']].sum() / (2 * ss['x']) - intertempelasticitybkmm
                else:
                    diffinter = G['micro']['realx']['pdur'][0:6,6:T['micro']].sum() / (ss['x'] * 6) - intertempelasticity
                
                
                
                if parameter=='gamma':
                    gamma = gamma - 0.1 * diffmpc
                   

                if ss['x']>10**-6:
                    
                    diff = np.max(np.abs([diffmpc, diffinter])) 
                    # thetad = np.max([thetad + 0.005 * diffinter,0])
                    if parasetname=='lesselastic':
                        thetad = np.max([thetad + 0.001 * diffinter,0])
                    else:
                        thetad = np.max([thetad + 0.005 * diffinter,0])
                    # diff = np.max(np.abs([diffmpc, diffinter])) 
                    print(diff)
                else:
                    diff = np.abs(diffmpc)


                if diff<10**-6:
                    break
                
            assert kk<maxit-1, 'did not converge'
            
            print('Found parameters')
            print(gamma, thetad)
            G['GE'] = jac.get_G(block_list=block_list['GE'],
                                exogenous=exog['GE'],
                                unknowns=unknowns['GE'],
                                targets=targets['GE'],
                                T=T['GE'], ss=ss)
            
            print('Solved GE')
            print(G['micro']['realx']['pdur'][0:6,6:T['micro']].sum() / (ss['x'] * 6))
            print(G['micro']['realx']['pdur'][0:12,12:T['micro']].sum() / (ss['x'] * 12))
            print(G['micro']['realx']['pdur'][0:lrhorizon,0:T['micro']].sum() / (ss['x'] * lrhorizon))
            
            

            # calculating MPCs (combine with later?)
            GEdict = {}
            for key in T.keys():
                vardict = {}
                for var in model_names['Real']:
                    hdict = {}
                    for h in range(12):
                        hdict[str(h+1) + ' month'] = G[key][var]['tshock'][h,0].sum() / ss['exp']
                    vardict[var] = hdict
                    hdict['3 months cumulative'] = G[key][var]['tshock'][0:3,0].sum() / ss['exp']
                    hdict['12 months cumulative'] = G[key][var]['tshock'][0:12,0].sum()  / ss['exp']                           
                GEdict[key] = vardict

            print('Computed MPC')

            # calculating IRFs
            modelshock = dict()

            for info in info_set: # 
                modelshockGE = dict()

                for GEkey in G.keys():
                    modelshockvar = dict()

                    for variable in G['GE'].keys():
                        
                    
                        if variable not in G[GEkey].keys():
                            continue

                        if info in ['full-info']:
                            modelshockbase = G[GEkey][variable]['tshock'][0:Tmax,0:Tmax] @ taxshock 
                        elif info in ['no-anticipation']:
                            modelshockbase = G[GEkey][variable]['tshock'][0:Tmax,0] * taxshock[0]
                            for t in range(1,Tmax):
                                modelshockbase[t:] = modelshockbase[t:] + G[GEkey][variable]['tshock'][0:Tmax-t,0] * taxshock[t]

                        modelshockvar[variable] = modelshockbase
                    
                    modelshockGE[GEkey] = modelshockvar
                    
                modelshock[info] = modelshockGE

                

            for info in modelshock.keys():
                for GE in modelshock[info].keys():
                    for variable in modelshock[info][GE].keys():
                    
                        if variable not in ratelist:
                            if ss[variable]!=0:
                                modelshock[info][GE][variable] = modelshock[info][GE][variable] / ss[variable]
                            else:
                                modelshock[info][GE][variable] = modelshock[info][GE][variable] / ss['y']
                        elif variable in ratelist:
                            modelshock[info][GE][variable] = (1 + modelshock[info][GE][variable])**12  - 1            
            print('Computed IRFs')

            mpc['micro-MPC = ' + str(round(mpctarget,2))] = {'MPCs': GEdict,
                                    'value': ss[parameter],
                                    'parameter': parameter,
                                    'steady': ss,
                                    'G': G,
                                    'shock': modelshock} 

        output[parameter] = mpc

    # save output
    pickle.dump( output , open( '../output/' + parasetname + '.pkl', 'wb' ) )

# %%
