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

from utilityfunction import f_u, f_invuc, f_uc, f_ud

from scipy.interpolate import LinearNDInterpolator, interpn
#%%

# TODO


def household_init(a_grid, d1_grid, netinc, r, p, eis, psi, deltad, f, vainit = 0, vdinit = 0, vinit = 0):
    coh = (1 + r) * a_grid[np.newaxis, np.newaxis, :] + p * (1 - deltad) * d1_grid[np.newaxis, :, np.newaxis] + netinc[:, np.newaxis, np.newaxis]
    Va = (1 + r) * (0.1 * psi * coh) ** (-1 / eis)
    Vd = (1-deltad) * (1 - f) * Va
    
    V = np.cumsum(Va, axis=2) + np.cumsum(Vd, axis=1)

    if vainit!=0:
        Va = vainit
        Vd = vdinit
        V = vinit
    
    return Va, Vd, V



@het(exogenous='Pi', policy=['d1','a'], backward=['V','Vd','Va'], backward_init=household_init)
def household_v(Va_p, Vd_p, V_p,  a_grid, d1_grid, e_grid, x_grid, m_grid, netinc, r, spread, p, user_cost, beta, eis, psi, deltad, xi, probadj, op_cost, maint_cost, collateral, d1_mesh, a_mesh, e_mesh, e_mesh_fine, x_mesh_fine, s_mesh_fine, nE, nSfine, nA, nX, nD, f, dmax, sigmaV, pbar, c_min):

    nN = np.prod(Va_p.shape)

    deltad_maint = deltad * (1 - maint_cost)

    maint_cost_ratio = maint_cost * deltad / (1 - deltad * (1 - maint_cost))

    noadj_dur_upkeep = pbar * op_cost + maint_cost_ratio - collateral * pbar

    effective_dur_price = p + pbar * op_cost - collateral * pbar

    # effective interest rate based on whether household was borrowing
    atilde_mesh = a_mesh - collateral * pbar * d1_mesh

    aborrow = np.zeros(a_mesh.shape)
    if  spread>0:
        borrowed = (atilde_mesh<0)
        interestrate = r * np.ones(atilde_mesh.shape)
        interestrate[borrowed] = r + spread
        aborrow[borrowed] = - atilde_mesh[borrowed]
    else:
        interestrate = r
        

    # marginal effect of existing durable stock on COH if not adjusting
    noadj_dur_value = noadj_dur_upkeep + collateral * pbar * (1 + interestrate) / ((1 - deltad_maint))

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
                  cohgrid = m_grid[np.newaxis, np.newaxis, :], 
                  dgrid = d1_grid[np.newaxis, :, np.newaxis],
                  dur_upkeep = noadj_dur_upkeep,
                  util_params = util_params):

        c_seq = f_invuc(beta * Va_p, dgrid, **util_params)

        coh0 = (savegrid 
            + c_seq 
            + dur_upkeep * dgrid)

        aegm = interpolate_y(coh0, cohgrid, savegrid)

        coh_net_cost = cohgrid - dur_upkeep * dgrid
        aegm = np.minimum(aegm, coh_net_cost)
        aegm = np.maximum(aegm, a_grid.min())

        return aegm

    aegm = solve_egm()

    # === STEP 2: Find no adjustment solution ===
    # no adjustment sets d'=(1-delta)d and has cash on hand x=a+y
    d1_noadj = (1 - deltad_maint) * d1_mesh
    coh_noadj = netinc[:, np.newaxis, np.newaxis] + (1 + interestrate) * atilde_mesh

    def solve_noadj(aegm = aegm, 
                    dprime = d1_noadj,
                    coh = coh_noadj, 
                    cohgrid = m_grid, 
                    dgrid = d1_grid,
                    agrid = a_grid, 
                    dur_upkeep = noadj_dur_upkeep,
                    cmin = c_min):

        i_d, pi_d = interpolate.interpolate_coord(dgrid, (1 - deltad_maint) * dgrid)
        i_m, pi_m = interpolate.interpolate_coord(cohgrid, coh)
        
        aprime = interpolate.apply_coord(i_d, pi_d, aegm.swapaxes(1, 2)).swapaxes(1, 2)
        aprime = interpolate.apply_coord(i_m, pi_m, aprime)

        coh_net_cost = coh - dur_upkeep * dprime
        aprime = np.minimum(aprime, coh_net_cost)
        aprime = np.maximum(aprime, agrid.min())
        
        cprime = ( coh_net_cost - aprime )

        # with borrowing there may not be a way to maintain positive consumption if not adjusting
        # in that case set consumption to a very small positive value
        cprime = np.maximum(cprime, cmin)

        return cprime, aprime

    c_noadj, a_noadj = solve_noadj()

     
    # === STEP 3: Update value function for no adjustment ===
    # V_a^noadj
    Va_noadj = (1 + interestrate) * f_uc(c_noadj, d1_noadj, **util_params)

    # V_d^noadj
    i_a, pi_a = interpolate.interpolate_coord(a_grid, a_noadj)
    i_d, pi_d = interpolate.interpolate_coord(d1_grid, (1 - deltad_maint) * d1_grid)
    
    Vd_p_inv = Vd_p ** -1
    Vd_p_inv_noadj = interpolate.apply_coord(i_d, pi_d, Vd_p_inv.swapaxes(1, 2)).swapaxes(1, 2)
    Vd_p_inv_noadj = interpolate.apply_coord(i_a, pi_a, Vd_p_inv_noadj)
    Vd_p_noadj = Vd_p_inv_noadj ** -1

    Vd_noadj =  (1 - deltad_maint) * (f_ud(c_noadj, d1_noadj, **util_params) 
                                - noadj_dur_value / (1 + interestrate) * Va_noadj  
                                + beta * Vd_p_noadj)

    V_p_noadj = interpolate.apply_coord(i_d, pi_d, V_p.swapaxes(1, 2)).swapaxes(1, 2)
    V_p_noadj = interpolate.apply_coord(i_a, pi_a, V_p_noadj)
    V_noadj = f_u(c_noadj, d1_noadj, **util_params) + beta * V_p_noadj

    # ============================
    # ===   ADJUSTER PROBLEM   ===
    # ============================  
    def solve_adjust(aegm = aegm,
                    coh_noy = x_mesh_fine,
                    income = netinc[:, np.newaxis, np.newaxis],
                    durexp_share = s_mesh_fine,
                    productivity = e_mesh_fine,
                    egrid = e_grid,
                    dgrid = d1_grid,
                    mgrid = m_grid,
                    agrid = a_grid,
                    interestrate = interestrate,
                    dur_value = adj_dur_value,
                    dur_price = effective_dur_price,
                    dur_upkeep = noadj_dur_upkeep,
                    dmax = dmax,
                    util_params = util_params):

        nNfine = np.prod(coh_noy.shape)
        nE,nX,nS = coh_noy.shape

        # === STEP 4: Create candidate values for optimal d' if adjust ===
        coh = coh_noy + income
        dchoice = np.minimum(coh / dur_price, dmax) * durexp_share
        mchoice = coh - (dur_price - dur_upkeep) * dchoice 

        interp_points = np.concatenate((productivity.reshape([nNfine,1]), dchoice.reshape([nNfine,1]), mchoice.reshape([nNfine,1])), axis=1)
        achoice = interpn((egrid, dgrid, mgrid), aegm, interp_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nX,nSfine])

        c_max = (mchoice - dur_upkeep * dchoice) 

        achoice = np.minimum(achoice, c_max - 10**-4)
        achoice = np.maximum(achoice, agrid.min())

        cchoice = c_max - achoice 

        # === STEP 5: Evaluate FOC for d' at candidate points ===
        interp_points[:,-1] = achoice.reshape([nNfine,])
        # candidate_points = np.concatenate((income.reshape([nNfine,1]), dchoice.reshape([nNfine,1]), achoice.reshape([nNfine,1])), axis=1)
        # assert np.max(np.abs(interp_points - candidate_points))<10**-8, 'interp'

        Vd_p_inv_interp = interpn((egrid, dgrid, agrid), Vd_p ** -1, interp_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nX,nSfine])
        Vd_p_interp = Vd_p_inv_interp ** -1

        evalfoc = ( - f_ud(cchoice, dchoice, **util_params)  
                    + dur_price * f_uc(cchoice, dchoice, **util_params)
                    - beta * Vd_p_interp )
        # print(evalfoc.shape)
        # === STEP 5b: FIND LOCATIONS WHERE FOC CHANGES SIGN ===
        focsign = np.concatenate((np.sign(evalfoc), np.ones([nE,nX,1])), axis=2)
        changesign = np.abs(np.diff(focsign, axis=2))
        print(changesign.sum(axis=2).max())
        # checks for multiple sign changes
        if changesign.sum(axis=2).max() > 2:
            # print(changesignloc)
            changesignsum = changesign.sum(axis=2)
            test = np.where(changesignsum>2)
            print(changesign[test])
            print(test)
            # WITH MULTIPLE INDECES NEED TO COMPARE OBJECTIVE FUNCTION
            # V_candidates = f_u(c_candidates, d1_candidates,  **util_params) + beta * V_p_candidates
            # s_index_max = np.argmax(V_candidates, axis=2, keepdims=True)
            STOP

        changesignloc = np.where(changesign)
        # STOP 

        changesignloc = np.minimum(changesignloc[2].reshape([nE,nX,1]), nS - 2)
        changesigncomb = np.concatenate((changesignloc,changesignloc+1), axis=2)

        changefoc = np.take_along_axis(evalfoc, changesigncomb, axis=2)
        changed = np.take_along_axis(dchoice, changesigncomb, axis=2)
        changec = np.take_along_axis(cchoice, changesigncomb, axis=2)


        # === STEP 6: Find solution to FOC by finding crossing point with 0 ===
        i_foc, pi_foc = interpolate.interpolate_coord(changefoc, np.zeros([nE,nX,1]))
        d_xadj = np.squeeze(interpolate.apply_coord(i_foc, pi_foc, changed))
        c_xadj = np.squeeze(interpolate.apply_coord(i_foc, pi_foc, changec))

        # limit solution to lie within interpolation bounds
        d_xadj = np.maximum(np.minimum(d_xadj, dchoice.max(axis=2)), changed.min(axis=2))
        c_xadj = np.maximum(np.minimum(c_xadj, cchoice.max(axis=2)), changec.min(axis=2))
        
        print(d_xadj.min())
        # === STEP 7: Interpolate decision rule and value function onto original a,d grid ===
        # fixed cost applies here
        coh_noy = dur_value * dgrid[np.newaxis,:,np.newaxis] + (1 + interestrate) * agrid[np.newaxis,np.newaxis,:]
        i_adj, pi_adj = interpolate.interpolate_coord(x_grid[np.newaxis,np.newaxis,:], coh_noy)
        d_adj = interpolate.apply_coord(i_adj, pi_adj, d_xadj[:,np.newaxis,:])
        c_adj = interpolate.apply_coord(i_adj, pi_adj, c_xadj[:,np.newaxis,:])
        
        print(d_adj.min())
        # coh = coh_noy + y[:, np.newaxis, np.newaxis]
        a_adj = (coh_noy + income - c_adj - dur_price * d_adj)
        a_adj = np.maximum(a_adj, agrid.min())
        
        return d_adj, a_adj, c_adj


    d1_adj, a_adj, c_adj = solve_adjust()

    # === STEP 8: Update value function for adjustment ===
    Va_adj = (1 + interestrate) * f_uc(c_adj, d1_adj, **util_params)
    Vd_adj = adj_dur_value / (1 + interestrate) * Va_adj 

    adjust_points = np.concatenate((e_mesh.reshape([nN,1]), d1_adj.reshape([nN,1]), a_adj.reshape([nN,1])), axis=1)
    V_p_adj = interpn((e_grid, d1_grid, a_grid), V_p, adjust_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nD,nA])
    V_adj = f_u(c_adj, d1_adj, **util_params) + beta * V_p_adj


    # === STEP 9: Adjustment function  ===
    meanshock = np.log(probadj / (1 - probadj))
    minprob = 30
    V_diff_scale = np.maximum(np.minimum((V_adj - V_noadj) / sigmaV + meanshock, minprob), -minprob)
    logadjust = V_diff_scale - np.log(1 + np.exp(V_diff_scale))
    adjust = np.exp(logadjust)

    # === STEP 10: Update solutions and value function ===
    a = a_adj * adjust + a_noadj * (1-adjust)
    c = c_adj * adjust + c_noadj * (1-adjust)
    d1 = d1_adj * adjust + d1_noadj * (1-adjust)
    Va = Va_adj * adjust + Va_noadj * (1-adjust)
    Vd = Vd_adj * adjust + Vd_noadj * (1-adjust)
    V = V_adj - sigmaV * logadjust

    # print(c_adj.min(), c_noadj.min())
    # print(d1_adj.min(), d1_noadj.min())

    anet = (a - collateral * pbar * d1 ) 
    
    assert np.isnan(V).sum()==0

    # check budget constraint
    coh_adj = (1 + interestrate) * a_mesh + adj_dur_value * d1_mesh + netinc[:,np.newaxis,np.newaxis]
    coh_noadj = (1 + interestrate) * atilde_mesh + netinc[:,np.newaxis,np.newaxis]
    assert np.max(np.abs(a_adj + effective_dur_price * d1_adj + c_adj - coh_adj) * adjust) < 10**-8, 'budget_adj'
    assert np.max(np.abs(a_noadj  + noadj_dur_upkeep * d1_noadj + c_noadj - coh_noadj) * (1-adjust) * (c_noadj > c_min)) < 10**-8, 'budget_noadj'


    expen = a_adj + c + effective_dur_price * d1 - (1-deltad) * p * d1_grid[np.newaxis,:,np.newaxis]
    incom = (1 + interestrate) * (a_mesh - collateral * pbar * d1_mesh) + netinc[:,np.newaxis,np.newaxis]
    saving = a_mesh
    fixedcost = f * (1 - deltad) * d1_mesh * adjust
    muc = Va

    
    
    return Va, Vd, V, a, c, d1, aborrow, anet, incom, expen, saving, fixedcost, muc, adjust
    
    # , a_noadj, c_noadj, d1_noadj, a_adj, c_adj, d1_adj, adjust
    
    # , v_adj, v_noadj, vd_adj, vd_noadj, va_adj, va_noadj





'''Part 1: Blocks'''

@simple
def firm(Y, Z_ss):
    w = Z_ss 
    L = Y / Z_ss
    return w, L

@simple
def monetary_policy(r_ss, rshock):
    r = r_ss + rshock
    return r  

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

@simple
def fiscal(G_ss, B_ss, r, Y, tau_ss, tshock, phib):
    G = G_ss
    Tr = tshock
    B = B_ss 
    tau =  ( G - B  + (1 + r) * B(-1) + Tr ) / Y
    return B, G, Tr, tau

@simple
def fiscal_ss(G_ss, B_ss, tau_ss, r, Y):
    G = G_ss
    Tr = 0
    B = B_ss 
    tau = tau_ss

    return B, G, Tr, tau 

@simple
def usercost(p, r, deltad, op_cost, pbar, spread):
    user_cost = p + pbar * op_cost - p(+1) * (1 - deltad) / (1 + r(+1))
    user_cost_borrow = p + pbar * op_cost - p(+1) * (1 - deltad) / (1 + r(+1) + spread)
    p_p = p(+1)
    return user_cost, user_cost_borrow, p_p


@simple
def mkt_clearing(B, ANET, ABORROW, FIXEDCOST, Y, C, G, p, D1, deltad, op_cost, pbar, spread, collateral): #, D, p
    asset_mkt = ANET - B 
    X = (D1 - (1 - deltad) * D1(-1) + FIXEDCOST)
    borrow_cost = spread * ABORROW
    goods_mkt = Y - G - C - p * X - pbar * op_cost * D1 - borrow_cost
    
    return asset_mkt, X, goods_mkt, borrow_cost
    # , cap_acc
    # , goods_mkt





'''Part 2: Embed HA block'''

def make_grids(rhoY, sigmaY, nE, amax, nA, dmax, nD, dmin, xmax, nX, xmin, nSfine, phiD, nM):
    e_grid, _, Pi = utils.discretize.markov_rouwenhorst(rho=rhoY, sigma=sigmaY, N=nE)
    a_grid = utils.discretize.agrid(amax=amax, n=nA)
    d1_grid = utils.discretize.agrid(amax=dmax, n=nD, amin=dmin)
    s_grid = np.linspace(start=0.001, stop=0.999, num=nSfine)
    x_grid = utils.discretize.agrid(amax=xmax, n=nX, amin=xmin)
    m_grid = utils.discretize.agrid(amax=xmax, n=nM, amin=xmin)

    e_mesh, d1_mesh, a_mesh = np.meshgrid(e_grid, d1_grid, a_grid, indexing='ij')
    e_mesh_ex, x_mesh_ex = np.meshgrid(e_grid, x_grid,  indexing='ij')
    e_mesh_fine, x_mesh_fine, s_mesh_fine = np.meshgrid(e_grid, x_grid, s_grid, indexing='ij')
    
    return e_grid, Pi, a_grid, d1_grid, s_grid, x_grid, m_grid, e_mesh, d1_mesh, a_mesh, e_mesh_fine, s_mesh_fine, x_mesh_fine, e_mesh_ex, x_mesh_ex


def income(w, tau, L, Tr, e_grid):

    netinc = w * (1 - tau) * e_grid * L + Tr

    return netinc


'''Part 3: DAG'''

def dag():
    # Combine blocks
    household = household_v.add_hetinputs([income, make_grids])
    fc_model = create_model([household, firm, usercost, mkt_clearing, fiscal, durableprice, monetary_policy], name="FC")
    fc_model_ss = create_model([household, firm_ss, usercost, mkt_clearing, fiscal_ss], name="FC SS")

    r = 0.01
    amax = 200
    dmax = 12
    xmax = amax + dmax
    B_ss = 1
    G_ss = 0.2
    tau_ss = G_ss + r * B_ss
    calibration = { 'eis': 1, 'xi': 1, 
                    'rhoY': 0.966, 'sigmaY': 0.5, 'sigmaV': 0.1,
                    'r': r, 'rshock': 0, 'spread': 0,
                    'f': 0.05, 'psi': 0.85, 'deltad': 1 - 0.8**(1/4), 
                    'op_cost': 0, 'maint_cost': 0, 'collateral': 0, 
                    'probadj': 0.5,
                    'pbar': 1, 'pshock': 1, 'psid': 0,
                    'Y': 1.0, 'L': 1.0, 
                    'G_ss': G_ss, 'B_ss': B_ss, 'tau_ss': tau_ss, 'phib': 0.1, 'tshock': 0,
                    'nE': 2, 'nA': 100, 'amax': amax, 'nD': 50, 'dmax': dmax, 'dmin': 0.05, 'phiD': 1, 
                    'nX': 200, 'xmax': xmax, 'xmin': 0, 'nSfine': 40, 'nM': 200,
                    'c_min': 10 ** -4} 
    unknowns_ss = {'beta': (0.975/(1+r), 0.995/(1+r))}
    targets_ss = {'asset_mkt': 0.}
    ss = fc_model_ss.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='brentq')

    # Transitional dynamics
    exogenous = ['tshock','rshock','pshock']
    unknowns = ['Y','p']
    targets = ['asset_mkt','dur_supply']
    ss['X_ss'] = ss['X']
    ss['r_ss'] = ss['r']
    ss['Z_ss'] = ss['Z']

    return fc_model_ss, ss, fc_model, unknowns, targets, exogenous

if __name__=='__main__':
    t0 = time.time()

    fc_model_ss, ss, fc_model, unknowns, targets, exogenous = dag()

    T = 300
    G = fc_model.solve_jacobian(ss, unknowns, targets, exogenous, T=T)

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
    plt.xlabel(r'quarters')
    plt.show()

    dZ = np.zeros([T,])
    dZ[0] = 0.01
    # dZ = 0.01*ss['Z']*rhos**(np.arange(T)[:, np.newaxis]) # get T*5 matrix of dZ
    dr = G['Y']['pshock'] @ dZ

    plt.plot(dr[:50, ])
    plt.title(r'$Y$ response to 1% $p$ shock$')
    plt.ylabel(r'basis points deviation from ss')
    plt.xlabel(r'quarters')
    plt.show()

    dZ = np.zeros([T,])
    dZ[0] = 0.01
    # dZ = 0.01*ss['Z']*rhos**(np.arange(T)[:, np.newaxis]) # get T*5 matrix of dZ
    dr = G['Y']['tshock'] @ dZ

    plt.plot(dr[:50, ])
    plt.title(r'$Y$ response to 1% $T$ shock$')
    plt.ylabel(r'basis points deviation from ss')
    plt.xlabel(r'quarters')
    plt.show()

    # tests

    for mkt in ['goods_mkt', 'asset_mkt']:
        assert ss[mkt]<10**-6, mkt + ' does not clear in steady state'
        assert abs(G[mkt]['rshock']).max()<10**-6, mkt + ' does not clear dynamically'


# %%
