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
from sequence_jacobian import simple, het, create_model, estimation, interpolate, misc
from sequence_jacobian.utilities.interpolate import interpolate_y

from utilityfunction import f_u, f_invuc, f_uc, f_ud, f_invuc_dc, f_uc_dc

from scipy.interpolate import LinearNDInterpolator, interpn
#%%

# TODO
# specify shock so it does not change unconditional probability of adjustment
# Income tax: do this in income block

def household_init(a_grid, d1_grid, y, r, p, eis, psi, deltad, f, vainit = 0, vdinit = 0, vinit = 0):
    coh = (1 + r) * a_grid[np.newaxis, np.newaxis, :] + p * (1 - deltad) * d1_grid[np.newaxis, :, np.newaxis] + y[:, np.newaxis, np.newaxis]
    Va = (1 + r) * (0.1 * psi * coh) ** (-1 / eis)
    Vd = (1-deltad) * (1 - f) * Va
    
    V = np.cumsum(Va, axis=2) + np.cumsum(Vd, axis=1)

    if vainit!=0:
        Va = vainit
        Vd = vdinit
        V = vinit
    
    return Va, Vd, V



@het(exogenous='Pi', policy=['d1','a'], backward=['V','Vd','Va'], backward_init=household_init)
def household_v(Va_p, Vd_p, V_p,  a_grid, d1_grid, s_grid, x_grid, m_grid, y, r, spread, p, user_cost, beta, eis, psi, deltad, xi, probadj, op_cost, maint_cost, collateral, d1_mesh, a_mesh, y_mesh, y_mesh_fine, x_mesh_fine, s_mesh_fine, nE, nSfine, nA, nX, nD, d1_candidates, x_candidates, f, dmax, sigmaV, pbar):

    nN = np.prod(y_mesh.shape)

    # === STEP 0: Define parameters (TBD: move out) ===

    deltad_maint = deltad * (1 - maint_cost)

    maint_cost_ratio = p * maint_cost * deltad / (1 - deltad_maint)
    
    
    noadj_dur_upkeep = pbar * op_cost + maint_cost_ratio

    effective_dur_price = p + pbar * op_cost - collateral * pbar

    # effective interest rate based on whether household was borrowing
    atilde_mesh = a_mesh - collateral * pbar * d1_mesh

    if  spread>0:
        borrowed = (atilde_mesh<0)
        interestrate = r * np.ones(atilde_mesh.shape)
        interestrate[borrowed] = r + spread
    else:
        interestrate = r

    noadj_dur_borrowing = collateral * pbar * (interestrate + deltad_maint) / (1 - deltad_maint)

    # marginal effect of existing durable stock on COH if not adjusting
    noadj_dur_value = noadj_dur_upkeep + noadj_dur_borrowing

    # marginal effect of existing durable stock on COH if adjusting
    adj_dur_value = (p * (1 - deltad) * (1 - f) - (1 + interestrate) * collateral * pbar)   

    util_params = {'eis': eis, 'xi': xi, 'psi': psi}
    
    # ============================
    # === NON-ADJUSTER PROBLEM ===
    # ============================
    # === STEP 1: Find a for given d' using EGM ===
    def solve_egm(Va_p = Va_p, 
                  beta = beta,
                  savegrid = a_grid[np.newaxis, np.newaxis, :], 
                  mgrid = m_grid[np.newaxis, np.newaxis, :], 
                  dgrid = d1_grid[np.newaxis, :, np.newaxis],
                  dur_upkeep = noadj_dur_upkeep,
                  util_params = util_params):

        c_seq = f_invuc(beta * Va_p, dgrid, **util_params)

        # m0 could be negative, how to fix? also add lambda p d' ?
        m0 = (savegrid
            + c_seq 
            + dur_upkeep * dgrid)

        aegm = interpolate_y(m0, mgrid, savegrid)

        m_net_cost = mgrid - dur_upkeep * dgrid
        aegm = np.minimum(aegm, m_net_cost)
        aegm = np.maximum(aegm, a_grid.min())

        return aegm

    aegm = solve_egm()


    # === STEP 2: Find no adjustment solution ===
    # no adjustment sets d'=(1-delta)d and has cash on hand x=a+y

    d_noadj = (1 - deltad_maint) * d1_mesh

    coh_noadj = y_mesh + (1 + interestrate) * atilde_mesh - noadj_dur_borrowing * d_noadj

    def solve_noadj(aegm = aegm, 
                    dprime = d_noadj,
                    coh = coh_noadj, 
                    cohgrid = m_grid, 
                    dgrid = d1_grid,
                    agrid = a_grid,  
                    dur_upkeep = noadj_dur_upkeep):

        i_d, pi_d = interpolate.interpolate_coord(dgrid, (1 - deltad_maint) * dgrid)
        i_m, pi_m = interpolate.interpolate_coord(cohgrid, coh)
        
        aprime = interpolate.apply_coord(i_d, pi_d, aegm.swapaxes(1, 2)).swapaxes(1, 2)
        aprime = interpolate.apply_coord(i_m, pi_m, aprime)

        coh_net_cost = coh - dur_upkeep * dprime
        aprime = np.minimum(aprime, coh_net_cost )
        aprime = np.maximum(aprime, agrid.min())
        
        cprime = coh_net_cost - aprime 

        # with borrowing there may not be a way to maintain positive consumption if not adjusting
        # in that case set consumption to a very small positive value
        cprime = np.maximum(cprime, 10 ** -4)

        return cprime, aprime

    c_noadj, a_noadj = solve_noadj()

    # === STEP 3: Update value function for no adjustment ===    
    Va_noadj = (1 + interestrate) * f_uc(c_noadj, d_noadj, **util_params)

    i_a, pi_a = interpolate.interpolate_coord(a_grid, a_noadj)
    i_d, pi_d = interpolate.interpolate_coord(d1_grid, (1 - deltad_maint) * d1_grid)

    Vd_noadj_p_inv_interp = interpolate.apply_coord(i_d, pi_d, (Vd_p ** -1).swapaxes(1, 2)).swapaxes(1, 2)
    Vd_noadj_p_inv_interp = interpolate.apply_coord(i_a, pi_a, Vd_noadj_p_inv_interp)
    Vd_noadj_p_interp = Vd_noadj_p_inv_interp ** -1

    V_noadj_p_interp = interpolate.apply_coord(i_d, pi_d, V_p.swapaxes(1, 2)).swapaxes(1, 2)
    V_noadj_p_interp = interpolate.apply_coord(i_a, pi_a, V_noadj_p_interp)

    Vd_noadj =  (1 - deltad_maint) * (f_ud(c_noadj, d_noadj, **util_params) 
                                        - noadj_dur_value * Va_noadj / (1 + interestrate) 
                                        + beta * Vd_noadj_p_interp)
    
    V_noadj = f_u(c_noadj, d_noadj, **util_params) + beta * V_noadj_p_interp

    # ============================
    # ===   ADJUSTER PROBLEM   ===
    # ============================  
    def solve_adjust(aegm = aegm,
                    coh_noy = x_mesh_fine,
                    income = y_mesh_fine,
                    durexp_share = s_mesh_fine,
                    ygrid = y,
                    dgrid = d1_grid,
                    mgrid = m_grid,
                    agrid = a_grid,
                    interestrate = interestrate,
                    dur_value = adj_dur_value,
                    dur_price = effective_dur_price,
                    dur_upkeep = noadj_dur_upkeep,
                    dmax = dmax,
                    util_params = util_params):

        nNfine = np.prod(income.shape)
        nE,nX,nS = income.shape

        # === STEP 4: Create candidate values for optimal d' if adjust ===
        coh = coh_noy + income
        dchoice = np.minimum(coh / dur_price, dmax) * durexp_share
        mchoice = coh - (dur_price - dur_upkeep) * dchoice 

        interp_points = np.concatenate((income.reshape([nNfine,1]), dchoice.reshape([nNfine,1]), mchoice.reshape([nNfine,1])), axis=1)
        achoice = interpn((ygrid, dgrid, mgrid), aegm, interp_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nX,nSfine])

        c_max = (mchoice - dur_upkeep * dchoice) 

        achoice = np.minimum(achoice, c_max - 10**-4)
        achoice = np.maximum(achoice, agrid.min())

        cchoice = c_max - achoice

        # === STEP 5: Evaluate FOC for d' at candidate points ===
        interp_points[:,-1] = achoice.reshape([nNfine,])
        # candidate_points = np.concatenate((income.reshape([nNfine,1]), dchoice.reshape([nNfine,1]), achoice.reshape([nNfine,1])), axis=1)
        # assert np.max(np.abs(interp_points - candidate_points))<10**-8, 'interp'

        Vd_p_inv_interp = interpn((ygrid, dgrid, agrid), Vd_p ** -1, interp_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nX,nSfine])
        Vd_p_interp = Vd_p_inv_interp ** -1

        evalfoc = ( - f_ud(cchoice, dchoice, **util_params)  
                    + dur_price * f_uc(cchoice, dchoice, **util_params)
                    - beta * Vd_p_interp )

        # === STEP 6: Find solution to FOC by finding crossing point with 0 ===
        i_foc, pi_foc = interpolate.interpolate_coord(evalfoc, np.zeros([nE,nX,1]))
        d_xadj = np.squeeze(interpolate.apply_coord(i_foc, pi_foc, dchoice))
        c_xadj = np.squeeze(interpolate.apply_coord(i_foc, pi_foc, cchoice))

        # limit solution to lie within interpolation bounds
        d_xadj = np.maximum(np.minimum(d_xadj, dchoice.max(axis=2)), dchoice.min(axis=2))
        c_xadj = np.maximum(np.minimum(c_xadj, cchoice.max(axis=2)), cchoice.min(axis=2))

        # === STEP 7: Interpolate decision rule and value function onto original a,d grid ===
        # fixed cost applies here
        coh_noy = dur_value * dgrid[np.newaxis,:,np.newaxis] + (1 + interestrate) * agrid[np.newaxis,np.newaxis,:]
        i_adj, pi_adj = interpolate.interpolate_coord(x_grid[np.newaxis,np.newaxis,:], coh_noy)
        d_adj = interpolate.apply_coord(i_adj, pi_adj, d_xadj[:,np.newaxis,:])
        c_adj = interpolate.apply_coord(i_adj, pi_adj, c_xadj[:,np.newaxis,:])
        

        # coh = coh_noy + y[:, np.newaxis, np.newaxis]
        a_adj = coh_noy + ygrid[:, np.newaxis, np.newaxis] - c_adj - dur_price * d_adj
        a_adj = np.maximum(a_adj, agrid.min())

        return d_adj, a_adj, c_adj

    d_adj, a_adj, c_adj = solve_adjust()

    # === STEP 8: Update value function for adjustment ===
    Va_adj = (1 + interestrate) * f_uc(c_adj, d_adj, **util_params)
    Vd_adj = adj_dur_value / (1 + interestrate) * Va_adj

    adjust_points = np.concatenate((y_mesh.reshape([nN,1]), d_adj.reshape([nN,1]), a_adj.reshape([nN,1])), axis=1)
    V_p_adj = interpn((y, d1_grid, a_grid), V_p, adjust_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nD,nA])
    V_adj = f_u(c_adj, d_adj, **util_params) + beta * V_p_adj


    # === STEP 9: Adjustment function  ===
    meanshock = np.log(probadj / (1 - probadj))
    minprob = 30
    V_diff_scale = np.maximum(np.minimum((V_adj - V_noadj) / sigmaV + meanshock, minprob), -minprob)
    logadjust = V_diff_scale - np.log(1 + np.exp(V_diff_scale))
    adjust = np.exp(logadjust)

    # === STEP 10: Update solutions and value function ===
    a = a_adj * adjust + a_noadj * (1-adjust)
    c = c_adj * adjust + c_noadj * (1-adjust)
    d1 = d_adj * adjust + d_noadj * (1-adjust)
    Va = Va_adj * adjust + Va_noadj * (1-adjust)
    Vd = Vd_adj * adjust + Vd_noadj * (1-adjust)
    V = V_adj - sigmaV * logadjust

    # Vd = Vd_adj
    # V = V_adj
    # Va = Va_adj
    # a = a_adj
    # c = c_adj
    # d1 = d1_adj
    # print(adjust.min(), adjust.max(), adjust.mean(), beta)
    # print(logadjust.min(), logadjust.max())
    # print(V_adj.min(), V_adj.max(),V_adj.mean(), np.isnan(V_adj).sum())
    # print(V_noadj.min(), V_noadj.max(),V_noadj.mean(), np.isnan(V_noadj).sum())
    # print(V.min(), V.max(),V.mean(), np.isnan(V).sum())
    # print(Va.min(), Va.max(), np.isnan(V).sum())
    # print(Va_adj.min(), Va_adj.max(), np.isnan(Va_adj).sum())
    # print(Va_noadj.min(), Va_noadj.max(), np.isnan(Va_noadj).sum())
    # print(Vd.min(), Vd.max(), Vd.mean(), np.isnan(Vd).sum())
    # print(Vd_adj.min(), Vd_adj.max(), Vd_adj.mean(), np.isnan(Vd_adj).sum())
    # print(Vd_noadj.min(), Vd_noadj.max(), Vd_noadj.mean(), np.isnan(Vd_noadj).sum())
    # print(d_noadj.min(), d_noadj.max(), np.isnan(d_noadj).sum())
    # print(d1_adj.min(), d1_adj.max(), d1_adj.mean(), beta)
    # print(d1_xadj.min(), d1_xadj.max(), d1_xadj.mean(), np.isnan(d1_xadj).sum())
    # print(d1_index.min(), d1_in/dex.max(), np.isnan(d1_index).sum())
    # print(d1.min(), d1.max(), d1.mean(), np.isnan(d1).sum())
    # print(c.min(), c.max(), np.isnan(c).sum())
    # print(c_noadj.min(), c_noadj.max(), np.isnan(c_noadj).sum())
    # print(c_adj.min(), c_adj.max(), np.isnan(c_adj).sum())
    # print(c_xadj.min(), c_xadj.max(), np.isnan(c_xadj).sum())
    # print(c_index.min(), c_index.max(), np.isnan(c_index).sum())
    # print(c_seq.min(), c_seq.max(), np.isnan(c_seq).sum())
    # print(cegm.min(), cegm.max(), np.isnan(cegm).sum())
    # print(aegm.min(), aegm.max(), np.isnan(aegm).sum())
    # print(c_seq.min(), c_seq.max(), np.isnan(c_seq).sum())
    # Vd_adj = Va_adj * (1 - deltad)
    # print(beta)


    assert np.isnan(V).sum()==0

    v_adj = V_adj
    v_noadj = V_noadj
    vd_adj = Vd_adj
    vd_noadj = Vd_noadj
    va_adj = Va_adj
    va_noadj = Va_noadj
    #
    
    return Va, Vd, V, a, c, d1, a_noadj, c_noadj, d_noadj, a_adj, c_adj, d_adj, adjust, v_adj, v_noadj, vd_adj, vd_noadj, va_adj, va_noadj





'''Part 1: Blocks'''

@simple
def firm(K, L, Z, alpha, delta, pbar, pshock):
    r = alpha * Z * (K(-1) / L) ** (alpha-1) - delta
    w = (1 - alpha) * Z * (K(-1) / L) ** alpha
    Y = Z * K(-1) ** alpha * L ** (1 - alpha)
    p = pbar * pshock
    return r, w, Y, p

@simple
def usercost(p, r, deltad, op_cost, pbar, spread):
    user_cost = p + pbar * op_cost - p(+1) * (1 - deltad) / (1 + r(+1))
    user_cost_borrow = p + pbar * op_cost - p(+1) * (1 - deltad) / (1 + r(+1) + spread)
    p_p = p(+1)
    r_p = r(+1)
    return user_cost, user_cost_borrow, p_p, r_p


@simple
def mkt_clearing(K, A, Y, C, p, D1, delta, deltad, op_cost, pbar, collateral): #, D, p
    asset_mkt = (A - collateral * pbar * D1 ) - K
    X = p * (D1 - (1 - deltad) * D1(-1))
    I = K - (1 - delta) * K(-1)
    goods_mkt = Y - C - X - I - pbar * op_cost * D1
    
    return asset_mkt, X, I, goods_mkt
    # , cap_acc
    # , goods_mkt


@simple
def firm_ss(r, Y, L, delta, alpha, deltad, pbar):
    '''Solve for (Z, K) given targets for (Y, r).'''
    rk = r + delta
    K = alpha * Y / rk
    Z = Y / K ** alpha / L ** (1 - alpha)
    w = (1 - alpha) * Z * (K / L) ** alpha
    p = pbar
    return K, Z, w, p


'''Part 2: Embed HA block'''

def make_grids(rhoY, sigmaY, nE, amax, nA, dmax, nD, dmin, xmax, nX, xmin, nSfine, phiD, nM):
    e_grid, _, Pi = utils.discretize.markov_rouwenhorst(rho=rhoY, sigma=sigmaY, N=nE)
    a_grid = utils.discretize.agrid(amax=amax, n=nA)
    d1_grid = utils.discretize.agrid(amax=dmax, n=nD, amin=dmin)
    # d_grid = utils.discretize.nonlinspace(amax=dmax, n=nD, phi=phiD, amin=dmin)
    # s_grid = utils.discretize.agrid(amax=0, n=101, amin=1)
    s_grid = np.linspace(start=0.001, stop=0.999, num=nSfine)
    x_grid = utils.discretize.agrid(amax=xmax, n=nX, amin=xmin)
    m_grid = utils.discretize.agrid(amax=xmax, n=nM, amin=xmin)
    
    return e_grid, Pi, a_grid, d1_grid, s_grid, x_grid, m_grid


def income(w, e_grid, d1_grid, a_grid, x_grid, s_grid):
    y = w * e_grid  # add income tax here
    y_mesh, d1_mesh, a_mesh = np.meshgrid(y, d1_grid, a_grid, indexing='ij')
    y_mesh_yx, x_mesh_yx = np.meshgrid(y, x_grid,  indexing='ij')
    y_mesh_fine, x_mesh_fine, s_mesh_fine = np.meshgrid(y, x_grid, s_grid, indexing='ij')

    d1_candidates = (x_grid[np.newaxis,:,np.newaxis] + y[:,np.newaxis,np.newaxis]  ) * s_grid[np.newaxis,np.newaxis,:] 
    x_candidates = (x_grid[np.newaxis,:,np.newaxis] + y[:,np.newaxis,np.newaxis]  ) * (1 - s_grid[np.newaxis,np.newaxis,:])

    return y, y_mesh, d1_mesh, a_mesh, y_mesh_fine, s_mesh_fine, x_mesh_fine, y_mesh_yx, x_mesh_yx, d1_candidates, x_candidates


'''Part 3: DAG'''

def dag():
    # Combine blocks
    household = household_v.add_hetinputs([income, make_grids])
    ks_model = create_model([household, firm, usercost, mkt_clearing], name="Krusell-Smith")
    ks_model_ss = create_model([household, firm_ss, usercost, mkt_clearing], name="Krusell-Smith SS")

    # Steady state
    r = 0.01
    amax = 200
    dmax = 12
    xmax = amax + dmax
    calibration = {'eis': 1, 'delta': 0.05, 'alpha': 0.11, 'rhoY': 0.966, 'sigmaY': 0.5, 'sigmaV': 0.075,
                   'psi': 0.87, 'deltad': 1 - 0.8**(1/4), 'xi': 1, 'probadj': 0.074,
                   'f': 0.15, 'op_cost': 0.055, 'maint_cost': 0.466, 'pbar': 1, 'r': r, 'spread': 0, 
                   'Y': 1.0, 'L': 1.0, 'nE': 2, 'nA': 50, 'amax': amax, 'nD': 30, 'dmax': dmax, 'dmin': 0.05, 'phiD': 1, 
                   'nX': 200, 'xmax': xmax, 'xmin': 0, 'nSfine': 40, 'nM': 200} 
    unknowns_ss = {'beta': (0.985/(1+r), 0.995/(1+r))}
    targets_ss = {'asset_mkt': 0.}
    ss = ks_model_ss.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='brentq')

    # Transitional dynamics
    exogenous = ['Z']
    unknowns = ['K']
    targets = ['asset_mkt']

    return ks_model_ss, ss, ks_model, unknowns, targets, exogenous

if __name__=='__name__':
    t0 = time.time()

    ks_model_ss, ss, ks_model, unknowns, targets, exogenous = dag()

    T = 200
    G = ks_model.solve_jacobian(ss, unknowns, targets, exogenous, T=T)

    t1 = time.time()
    print(t1-t0)

    rhos = np.array([0.2, 0.4])
    dZ = 0.01*ss['Z']*rhos**(np.arange(T)[:, np.newaxis]) # get T*5 matrix of dZ
    dr = G['r']['Z'] @ dZ

    plt.plot(10000*dr[:50, :])
    plt.title(r'$r$ response to 1% $Z$ shocks with $\rho=(0.2 ... 0.9)$')
    plt.ylabel(r'basis points deviation from ss')
    plt.xlabel(r'quarters')
    plt.show()

    # dZ = 0.01*(np.arange(T)[:, np.newaxis] == np.array([5, 10]))
    # dK = G['K']['Z'] @ dZ

    # plt.plot(dK[:50])
    # plt.title('$K$ response to 1% Z news shocks for $t=5,...,25$')
    # plt.show()

    # tests

    for mkt in ['goods_mkt', 'asset_mkt', 'dur_acc']:
        assert ss[mkt]<10**-6, mkt + ' does not clear in steady state'
        assert abs(G[mkt]['Z']).max()<10**-6, mkt + ' does not clear dynamically'


