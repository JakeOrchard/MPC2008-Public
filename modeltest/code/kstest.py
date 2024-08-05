
#%%
import sys
sys.path.append('sequence_jacobian')

import copy
import numpy as np
from numba import njit
import scipy.optimize as opt
import scipy.linalg as linalg
import matplotlib.pyplot as plt

from sequence_jacobian import simple, het, create_model, estimation
import sequence_jacobian.utilities as utils

# import kstesthh as hh


def household_init(a_grid, y, r, eis):
    coh = (1 + r) * a_grid[np.newaxis, :] + y[:, np.newaxis]
    Va = (1 + r) * (0.1 * coh) ** (-1 / eis)
    return Va


@het(exogenous='Pi', policy='a', backward='Va', backward_init=household_init)
def household_v(Va_p, a_grid, y, r, beta, eis):
    uc_nextgrid = beta * Va_p
    c_nextgrid = uc_nextgrid ** (-eis)
    coh = (1 + r) * a_grid[np.newaxis, :] + y[:, np.newaxis]
    a = utils.interpolate.interpolate_y(c_nextgrid + a_grid, coh, a_grid)
    utils.optimized_routines.setmin(a, a_grid[0])
    c = coh - a
    Va = (1 + r) * c ** (-1 / eis)
    return Va, a, c


'''Part 1: Blocks'''

@simple
def firm(K, L, Z, alpha, delta):
    r = alpha * Z * (K(-1) / L) ** (alpha-1) - delta
    w = (1 - alpha) * Z * (K(-1) / L) ** alpha
    Y = Z * K(-1) ** alpha * L ** (1 - alpha)
    return r, w, Y


@simple
def mkt_clearing(K, A, Y, C, delta):
    asset_mkt = A - K
    goods_mkt = Y - C - delta * K
    return asset_mkt, goods_mkt


@simple
def firm_ss(r, Y, L, delta, alpha):
    '''Solve for (Z, K) given targets for (Y, r).'''
    rk = r + delta
    K = alpha * Y / rk
    Z = Y / K ** alpha / L ** (1 - alpha)
    w = (1 - alpha) * Z * (K / L) ** alpha
    return K, Z, w


'''Part 2: Embed HA block'''

def make_grids(rho, sigma, nS, amax, nA):
    e_grid, _, Pi = utils.discretize.markov_rouwenhorst(rho=rho, sigma=sigma, N=nS)
    a_grid = utils.discretize.agrid(amax=amax, n=nA)
    return e_grid, Pi, a_grid


def income(w, e_grid):
    y = w * e_grid
    return y


'''Part 3: DAG'''

def dag():
    # Combine blocks
    household = household_v.add_hetinputs([income, make_grids])
    ks_model = create_model([household, firm, mkt_clearing], name="Krusell-Smith")
    ks_model_ss = create_model([household, firm_ss, mkt_clearing], name="Krusell-Smith SS")

    # Steady state
    calibration = {'eis': 1.0, 'delta': 0.025, 'alpha': 0.11, 'rho': 0.966, 'sigma': 0.5,
                   'Y': 1.0, 'L': 1.0, 'nS': 2, 'nA': 10, 'amax': 200, 'r': 0.01}
    unknowns_ss = {'beta': (0.98 / 1.01, 0.999 / 1.01)}
    targets_ss = {'asset_mkt': 0.}
    ss = ks_model_ss.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='brentq')

    # Transitional dynamics
    exogenous = ['Z']
    unknowns = ['K']
    targets = ['asset_mkt']

    return ks_model_ss, ss, ks_model, unknowns, targets, exogenous

ks_model_ss, ss, ks_model, unknowns, targets, exogenous = dag()

T = 200
G = ks_model.solve_jacobian(ss, unknowns, targets, exogenous, T=T)


rhos = np.array([0.2, 0.4, 0.6, 0.8, 0.9])
dZ = 0.01*ss['Z']*rhos**(np.arange(T)[:, np.newaxis]) # get T*5 matrix of dZ
dr = G['r']['Z'] @ dZ

plt.plot(10000*dr[:50, :])
plt.title(r'$r$ response to 1% $Z$ shocks with $\rho=(0.2 ... 0.9)$')
plt.ylabel(r'basis points deviation from ss')
plt.xlabel(r'quarters')
plt.show()

dZ = 0.01*(np.arange(T)[:, np.newaxis] == np.array([5, 10, 15, 20, 25]))
dK = G['K'] @ dZ

plt.plot(dK[:50])
plt.title('$K$ response to 1% Z news shocks for $t=5,...,25$')
plt.show()

# G = rbc_model.solve_jacobian(ss, unknowns, targets, exogenous, T=T)

# %%
