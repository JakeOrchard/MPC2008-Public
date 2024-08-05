#%%
from multiprocessing.resource_sharer import stop
import sys
sys.path.append('sequence_jacobian')

import copy
import numpy as np
# from numba import njit
import scipy.sparse as sp
import scipy.optimize as opt
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import time

import sequence_jacobian.utilities as utils
from newton_solver2 import newton_solver2
from sequence_jacobian import simple, solved, het, create_model, estimation, interpolate, misc
from sequence_jacobian.utilities.interpolate import interpolate_y

from utilityfunction import f_u, f_invuc, f_uc, f_ud, f_invuc_dc, f_uc_dc

from scipy.interpolate import LinearNDInterpolator, interpn
#%%
from basicfixedcost import household_v, income, make_grids, household_init


# TODO
# measure MPC
# target durable size
# measure fraction of adjsutments from fixed costs?

'''Part 1: Blocks'''

@simple
def firm(Y, Z_ss):
    w = Z_ss 
    L = Y / Z_ss
    return w, L

@solved(unknowns={'sw': (0.9, 1.1), 'f1': (1, 100), 'f2': (1, 100), 'wstar': (0.5, 2)},
        targets=["distortw", "focw1", "focw2", 'focw'])
def union(MUC, L, sw, f1, f2, w, wstar, pi, beta, nu, phi, thetaw, epsw, chiw):
    
    H = sw * L
    distortw = (-sw + (1 - thetaw) * (wstar / w)**(-epsw) 
                + thetaw * (w(-1 ) /w)**(-epsw) * (pi / pi(-1)**chiw)**epsw * sw(-1) )

    focw1 = (-f1 + (epsw - 1) / epsw * wstar * MUC * (w / wstar)**epsw * L 
            + thetaw * beta * (pi(+1) / pi**chiw)**(epsw - 1) * (wstar(+1) / wstar)**(epsw - 1) * f1(+1) )
    focw2 = (-f2 + nu * H**phi * (w / wstar)**epsw * L 
            + thetaw * beta * (pi(+1) / pi**chiw)**epsw * (wstar(+1) / wstar)**epsw * f2(+1) )
    focw = -f1 + f2

    wages = (-w**(1-epsw) + (1-thetaw) * wstar**(1-epsw) 
            + thetaw * w(-1)**(1-epsw) * (pi(-1)**chiw / pi)**(1 - epsw) )        

    return H, focw1, focw2, focw, distortw, wages 

@simple
def union_ss(MUC, L, w, beta, phi, thetaw, epsw):
    
    sw = 1
    H = sw * L
    wstar = w
    f1 = (epsw - 1) / epsw * wstar * MUC * L / (1 - thetaw * beta)
    f2 = f1
    nu = f2 * (1 - thetaw * beta) / ( H**phi * (w / wstar)**epsw * L  )
    pi = 1 

    return sw, H, wstar, f1, f2, nu, pi    

# 
# @simple
@solved(unknowns={'r': (0.9, 1.1)}, targets=['realrate'])
def monetary_policy(r, pi, Y, r_ss, Y_ss, rshock, rhor, phipi, phigap):

    # intrule = -nomrate + rhor * nomrate(-1) + (1 - rhor) * (r_ss + phipi * (pi - 1) + phigap * (Y / Y_ss - 1)) 

    # realrate = -r + (1 + nomrate) / pi(+1) - 1  + rshock

    # return realrate, intrule 

    realrate = -r + rhor * r(-1) + (1 - rhor) * (r_ss + rshock + phipi * (pi - 1) + phigap * (Y / Y_ss - 1))

    return realrate

@simple
def durableprice(p, pbar, pshock, X, X_ss, psid):
    dur_supply = p - pbar * pshock * (X / X_ss) ** psid
    return dur_supply   

@simple
def firm_ss(Y, L, pbar):
    '''Solve for (Z, K) given targets for (Y, r).'''
    Z = Y / L
    w = Z
    p = pbar
    return w, Z, p

@solved(unknowns={'B': (0.1, 2)}, targets=["govtbudget"])
def fiscal(G_ss, B_ss, r, Y, B, tau_ss, tshock, phib):
    G = G_ss
    Tr = tshock
    tau = tau_ss + phib * (B(-12) - B_ss) 
    govtbudget = -G + B - (1 + r) * B(-1) - Tr + tau * Y

    return govtbudget, G, Tr, tau

@simple
def fiscal_ss(G_ss, B_ss, tau_ss, r, Y):
    G = G_ss
    Tr = 0
    B = B_ss 
    tau = tau_ss
    govtbudget =  -G + B - (1 + r) * B(-1) - Tr + tau * Y

    return B, G, Tr, tau, govtbudget 

@simple
def usercost(p, r, deltad, op_cost, pbar, spread):
    user_cost = p + pbar * op_cost - p(+1) * (1 - deltad) / (1 + r(+1))
    user_cost_borrow = p + pbar * op_cost - p(+1) * (1 - deltad) / (1 + r(+1) + spread)
    p_p = p(+1)
    return user_cost, user_cost_borrow, p_p


@simple
def mkt_clearing(B, ANET, ABORROW, FIXEDCOST, Y, C, G, p, D1, deltad, op_cost, pbar, spread, xc_ss): #, D, p
    asset_mkt = ANET - B 
    X = (D1 - (1 - deltad) * D1(-1) + FIXEDCOST)
    borrow_cost = spread * ABORROW
    goods_mkt = Y - G - C - p * X - pbar * op_cost * D1 - borrow_cost
    dur_target = xc_ss - X / (C + op_cost * D1 + borrow_cost)
    
    return asset_mkt, X, goods_mkt, borrow_cost, dur_target
    # , cap_acc
    # , goods_mkt


'''Part 3: DAG'''

if __name__=='__main__':
    # Combine blocks
    household = household_v.add_hetinputs([income, make_grids])
    hank_model = create_model([household, firm, usercost, mkt_clearing, fiscal, durableprice, monetary_policy, union], name="FC")
    hank_model_MPC = create_model([household, durableprice, mkt_clearing], name="FC MPC")
    hank_model_ss = create_model([household, firm_ss, usercost, mkt_clearing, fiscal_ss, union_ss], name="FC SS")

    r = 0.99**(-1/3) - 1 
    # r = 0.01
    amax = 200
    dmax = 15
    xmax = amax + dmax
    B_ss = 1
    G_ss = 0.2
    tau_ss = G_ss + r * B_ss
    rhoYq = 0.966
    rhoYm = 0.966**(1/3)
    sigmaYq = 0.5
    sigmaYm = sigmaYq * 3 / np.sqrt(1 + (1 + rhoYq) ** 2 + (1 + rhoYq + rhoYq**2) ** 2 + rhoYq ** 2 * (1 + rhoYq) ** 2  + rhoYq ** 2)
    calibration = { 'eis': 0.5, 'xi': 0.5, 
                    'rhoY': rhoYm, 'sigmaY': sigmaYm, 'sigmaV': 0.25,
                    'r': r, 'rshock': 0, 
                    'f': 0.01, 'psi': 0.75, 'deltad': 1 - 0.8**(1/12), 
                    'op_cost': 0.055/3, 'maint_cost': 0.466,
                    'collateral': 0.8, 'spread': r,
                    'probadj': 0.025,
                    'pbar': 1, 'pshock': 1, 'psid': 0,
                    'phi': 1,                   # 1/phi is Frisch elasticity) GLV =  0.2
                    'epsw': 6,                  # wage elasticity
                    'thetaw': 0.9167,           # wage stickiness, monthly equiv to theta = 0.75 at quarterly, 0 implies no stickiness
                    'chiw': 0,                  # wage indexing, , 0=no, 1=full
                    'phipi': 0.5,               # Taylor rule parameter on inflation
                    'phigap': 1/12,             # Taylor rule parameter on output gap, 1 for annual in balanced-approach rule, divide by 12 for monthly
                    'rhor': 0.85**(1/3),        # Inertia in Taylor rule, 0 = no inertia, quarterly rhor = monthly rhor^3, .85 from Fed website
                    'xc_ss': 0.0411,            # ratio of durable to nondurable expenditure
                    'Y': 1.0, 
                    'L': 1.0, 
                    'G_ss': G_ss, 'B_ss': B_ss, 'tau_ss': tau_ss, 'phib': 0.1, 'tshock': 0,
                    'nE': 2, 'nA': 100, 'amax': amax, 'nD': 50, 'dmax': dmax, 'dmin': 0.05, 'phiD': 1, 
                    'nX': 200, 'xmax': xmax, 'xmin': 0, 'nSfine': 40, 'nM': 200,
                    'c_min': 10 ** -4} 
    unknowns_ss = {'beta': (0.9875/(1+r), 0.995/(1+r))}
    targets_ss = {'asset_mkt': 0.}

    t0 = time.time()
    ss = hank_model_ss.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='brentq')

    ss['X_ss'] = ss['X']
    ss['r_ss'] = ss['r']
    ss['Z_ss'] = ss['Z']
    ss['Y_ss'] = ss['Y']
    ss['nomrate'] = ss['r_ss']

    # Transitional dynamics
    T = 300
    # exogenous = ['tshock','rshock','pshock']
    # unknowns = ['Y','p','pi']
    # targets = ['asset_mkt','dur_supply','wages']
    
    # GMPC
    GMPC = hank_model_MPC.solve_jacobian(ss, 
                                  ['Y','p'], 
                                  ['goods_mkt','dur_supply'], 
                                  ['Tr','pshock'], 
                                  T=T)
    
    print(GMPC['Y']['Tr'][0:4,0:4])

    # GE
    G = hank_model.solve_jacobian(ss, 
                                  ['Y','p','pi'], 
                                  ['asset_mkt','dur_supply','wages'], 
                                  ['tshock','rshock','pshock'], 
                                  T=60)

    t1 = time.time()
    print(t1-t0)

    # rhos = np.array([0.2, 0.4])
    dZ = np.zeros([T,])
    dZ[1] = 0.01
    # dZ = 0.01*ss['Z']*rhos**(np.arange(T)[:, np.newaxis]) # get T*5 matrix of dZ
    dr = G['Y']['rshock'] @ dZ

    plt.plot(dr[:50, ])
    plt.title(r'$Y$ response to 1% $r$ shock$')
    plt.ylabel(r'basis points deviation from ss')
    plt.xlabel(r'months')
    plt.show()

    dZ = np.zeros([T,])
    dZ[0] = 0.01
    # dZ = 0.01*ss['Z']*rhos**(np.arange(T)[:, np.newaxis]) # get T*5 matrix of dZ
    dr = G['Y']['pshock'] @ dZ

    plt.plot(dr[:50, ])
    plt.title(r'$Y$ response to 1% $p$ shock$')
    plt.ylabel(r'basis points deviation from ss')
    plt.xlabel(r'months')
    plt.show()

    dZ = np.zeros([T,])
    dZ[0] = 0.01
    # dZ = 0.01*ss['Z']*rhos**(np.arange(T)[:, np.newaxis]) # get T*5 matrix of dZ
    dr = G['Y']['tshock'] @ dZ

    plt.plot(dr[:50, ])
    plt.title(r'$Y$ response to 1% $T$ shock$')
    plt.ylabel(r'basis points deviation from ss')
    plt.xlabel(r'months')
    plt.show()

    for mkt in ['goods_mkt', 'asset_mkt']:
        assert ss[mkt]<10**-6, mkt + ' does not clear in steady state'
        assert abs(G[mkt]['rshock']).max()<10**-6, mkt + ' does not clear dynamically'



# %%
