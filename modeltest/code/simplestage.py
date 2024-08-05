#%%

# TODO
# try simple durable choice where adjust = d / (1 - delta) and no adjust = d * (1 - delta)


import numpy as np
from numba import njit
import matplotlib.pyplot as plt
from copy import deepcopy
# import sequence_jacobian as sj

from sequence_jacobian import grids, interpolate
from sequence_jacobian.blocks.stage_block import StageBlock
from sequence_jacobian.blocks.support.stages import Continuous1D, Continuous2D, ExogenousMaker, LogitChoice

from scipy.interpolate import LinearNDInterpolator, interpn

prod_stage = ExogenousMaker(markov_name='z_markov', index=0, name='prod')



@njit
def util_c(c):
    u = np.log(c)
    return u

@njit
def margutil_c(c):
    muc = c ** -1
    return muc  

@njit
def margutilinv_c(muc):
    c = muc ** -1
    return c            

@njit
def upperenv_vec(W, m_endo, coh, m_grid, a_grid):
    
    """Interpolate value function and consumption to exogenous grid."""
    n_b, n_m = coh.shape
    n_b2, n_a = W.shape

    assert coh.min()>0

    assert np.abs(n_b-n_b2)<10**-8

    a = np.zeros_like(coh)
    c = np.zeros_like(coh)
    V = -np.inf * np.ones_like(coh)
    Vm = -np.inf * np.ones_like(coh)

    # loop over other states, collapsed into single axis
    for ib in range(n_b):
        # loop over segments of endogenous asset grid from EGM (not necessarily increasing)
        for ja in range(n_a - 1):
            m_low, m_high = m_endo[ib, ja], m_endo[ib, ja + 1]
            W_low, W_high = W[ib, ja], W[ib, ja + 1]
            ap_low, ap_high = a_grid[ja], a_grid[ja + 1]
            
           # loop over exogenous asset grid (increasing) 
            for im in range(n_m):  
                mcur = m_grid[im]
                coh_cur = coh[ib, im]
                
                interp = (m_low <= mcur <= m_high) 
                extrap = (ja == n_a - 2) and (mcur > m_endo[ib, n_a - 1])

                # exploit that a_grid is increasing
                if (m_high < mcur < m_endo[ib, n_m - 1]):
                    break

                if interp or extrap:
                    # is this the right extrap for W0, Wd0?
                    # W is defined over a grid, here we use m-grid to extrap...
                    W0 = interpolate.interpolate_point(mcur, m_low, m_high, W_low, W_high)
                    a0 = interpolate.interpolate_point(mcur, m_low, m_high, ap_low, ap_high)
                    c0 = coh_cur - a0
                    V0 = util_c(c0) + W0
                    Vm0 = margutil_c(c0)

                    # upper envelope, update if new is better
                    if V0 > V[ib, im]:
                        a[ib, im] = a0 
                        c[ib, im] = c0
                        V[ib, im] = V0
                        Vm[ib, im] = Vm0

        # Enforce borrowing constraint
        im = 0
        while im < n_m and m_grid[im] <= m_endo[ib, 0]:
            a[ib, im] = a_grid[0]
            c[ib, im] = coh[ib, im]
            V[ib, im] = util_c(c[ib, im]) + W[ib, 0]
            Vm[ib, im] = margutil_c(c[ib, im])
            im += 1

    

    return V, Vm, c, a                           

def upperenv(W, m_endo, coh, m_grid, a_grid):
    # collapse (n, z, a) into (b, a)
    shapea = W.shape
    W = W.reshape((-1, shapea[-1]))
    m_endo = m_endo.reshape((-1, shapea[-1]))

    shapem = coh.shape
    coh = coh.reshape((-1, shapem[-1]))

    V, Vm, c, a = upperenv_vec(W, m_endo, coh, m_grid, a_grid)

    # report on (n, z, a)
    return V.reshape(shapem), Vm.reshape(shapem), c.reshape(shapem), a.reshape(shapem)                           


def dcegm_func(uc_endo, a_grid, coh_m):
    """DC-EGM algorithm"""
    c_endo = margutilinv_c(uc_endo)                # Euler equation
    
    aegm = interpolate.interpolate_y(c_endo + a_grid[np.newaxis,:], coh_m, a_grid)
    aegm = np.maximum(aegm, a_grid[0])

    cegm = coh_m - aegm

    return aegm, cegm   

def dcegm(V, Va, a_grid, m_grid, r, y, coh_m, beta):
    """DC-EGM algorithm"""
    # use all FOCs on endogenous grid
    W = beta * V                                                  # end-of-stage vfun
    uc_endo = beta * Va                                           # envelope condition
    c_endo = margutilinv_c(uc_endo) 
    m_endo = c_endo + a_grid[np.newaxis,:] - y
    
    
    aegm, cegm = dcegm_func(uc_endo, a_grid, coh_m)

    Vegm, Vmegm, cegm2, aegm2 = upperenv(W, m_endo, coh_m, m_grid, a_grid)
    # print(aegm[0,:5])
    # print(aegm2[0,:5])
    assert np.max(np.abs(aegm - aegm2))<10 ** -8

    coh = (1 + r) * a_grid[np.newaxis, :] + y
    a = interpolate.interpolate_y(coh_m, coh, aegm)
    a = np.maximum(a, a_grid[0])


    c = coh - a
    Va = (1 + r) * margutil_c(c) 

    Vatest = (1 + r) * interpolate.interpolate_y(coh_m, coh, Vmegm)
    print(np.max(np.abs((Vatest - Va)/Va)))

    Vprime = interpolate.interpolate_y(a_grid, a, W)
    V = util_c(c) + Vprime

    Vtest = interpolate.interpolate_y(coh_m, coh, Vegm)
    print(np.max(np.abs(Vtest - V)))

    # print(V[0,-5:],Vtest[0,-5:],Vegm[0,-5:])

    # Va = Vatest
    # V = Vtest

    return V, Va, a, c   

consav_stage = Continuous1D(backward=['Va','V'], policy='a', f=dcegm,
                            name='consav')     



def hh_init(coh, a_grid):
    V = util_c(0.2 * coh) / 0.01    
    Va = np.empty_like(V)
    
    Va[..., 1:-1] = (V[..., 2:] - V[..., :-2]) / (a_grid[2:] - a_grid[:-2])
    Va[..., 0] = (V[..., 1] - V[..., 0]) / (a_grid[1] - a_grid[0])
    Va[..., -1] = (V[..., -1] - V[..., -2]) / (a_grid[-1] - a_grid[-2])
    # Vd = np.zeros_like(V)


    
    return V, Va
    #, Vd

def make_grids(rho_z, sd_z, n_z, min_a, max_a, n_a, min_m, max_m, n_m):
    z_grid, z_dist, z_markov = grids.markov_rouwenhorst(rho_z, sd_z, n_z)
    a_grid = grids.agrid(max_a, n_a, min_a)
    
    m_grid = grids.agrid(max_m, n_m, min_m)

    aegm_grid = m_grid
    return z_grid, z_dist, z_markov, a_grid, m_grid, aegm_grid


def labor_income(m_grid, a_grid, z_grid, r, w):
    y_grid = z_grid * w 

    y = y_grid[:, np.newaxis]           # on (n, z)
    coh_m = m_grid[np.newaxis, :] + y 
    coh = (1 + r) * a_grid[np.newaxis, :] + y 
    

    return y, coh_m, coh

hh = dict()
hh['joint'] = StageBlock([prod_stage, consav_stage], name='hh',
                backward_init=hh_init, hetinputs=[make_grids, labor_income])

               

# print(hh)
# print(f"Inputs: {hh.inputs}")
# print(f"Outputs: {hh.outputs}")

cali = dict()

cali = {'taste_shock': 1E-10, 'r': 0.02/4, 'beta': 0.97, 'eis': 1, 
               'psi': 0, 'deltad': 0.5,
               'rho_z': 0.95, 'sd_z': 0.5, 'n_z': 2, 
               'min_a': 0.0, 'max_a': 200, 'n_a': 100,
               'min_m': 0.0, 'max_m': 200, 'n_m': 400,
               'min_d': 0.01, 'max_d': 0.1, 'n_d': 10,
               'w': 1.0}

ss = dict()
for i in hh.keys():
    ss[i] = hh[i].steady_state(cali)
    print(f"Aggregate assets in model {i}: {ss[i]['A']:0.2f}")

def fig1(ss, amax=150, amin=0, iz=0, figsize=0.6):
    a_grid = ss['joint'].internals['hh']['a_grid']
    a, da, c, P, V = dict(), dict(), dict(), dict(), dict()
    models = ['joint']
    for i in models:
        a[i] = ss[i].internals['hh']['consav']['a']
        da[i] = a[i] - a_grid
        c[i] = ss[i].internals['hh']['consav']['c']
        V[i] = ss[i].internals['hh']['consav']['V']

    fig, axes = plt.subplots(1, 2, figsize=(12*figsize, 4*figsize))
    ax = axes.flatten()
    
    labels = [r'small taste shock', r'large taste shock'] 

    for i, l in zip(models, labels):
        ax[0].plot(a_grid[:amax], a[i][iz, :amax], label=l, linewidth=2)
    ax[0].plot(a_grid[:amax], a_grid[:amax], color='gray', linestyle=':')
    ax[0].legend(frameon=False)
    ax[0].axhline(0, color='gray', linestyle=':')
    ax[0].set_title('Assets')

    for i, l in zip(models, labels):
        ax[1].plot(a_grid[:amax], c[i][iz, :amax], label=l, linewidth=2)
    ax[1].set_title('Consumption')
    
    for k in ax:
        k.set_xlabel('assets')

    plt.tight_layout()
    plt.show()

fig1(ss, amax=85, iz=0, figsize=0.8)    
# %%
