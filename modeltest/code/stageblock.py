#%%


import numpy as np
from numba import njit
import matplotlib.pyplot as plt
from copy import deepcopy
# import sequence_jacobian as sj

from sequence_jacobian import grids, interpolate, simple, create_model
from sequence_jacobian.blocks.stage_block import StageBlock
from sequence_jacobian.blocks.support.stages import Continuous1D, Continuous2D, ExogenousMaker, LogitChoice
from sequence_jacobian.blocks.support.law_of_motion import DiscreteChoice
from sequence_jacobian.utilities.misc import logit_choice

from scipy.interpolate import LinearNDInterpolator, interpn

from func_upper_envelope import upperenv, upperenv_durable
from func_utility import util, margutil_c, margutil_d, margutilinv_c

prod_stage = ExogenousMaker(markov_name='z_markov', index=1, name='prod')


adjust_stage = LogitChoice(value='V', backward=['Va', 'Vd'], index=0, name='adjust',
                           taste_shock_scale='taste_shock')


def optimal_ad(V, Va, Vd, a_grid, m_grid, d1_grid, y_grid, x_grid, coh, interestrate, deltad_maint, dur_value, dur_price, beta, psi, eis, eisd):

    # parameters for utility function
    util_para = {'psi': psi, 'eis': eis, 'eisd': eisd}

    assert np.allclose(V[1,...],V[0,...])
    assert np.allclose(Vd[1,...],Vd[0,...])
    assert np.allclose(Va[1,...],Va[0,...])

    # lagged adjustment indicator irrelevant given current a,d,y, so drop
    # first dimension
    W = beta * V[1,...].squeeze()                    # end-of-stage vfun
    Wd = beta * Vd[1,...].squeeze()                  # end-of-stage vfun
    uc_endo = beta * Va[1,...].squeeze()             # envelope condition

    d1egm = d1_grid[np.newaxis,:,np.newaxis] + np.zeros_like(W)

    c_egm_agrid = margutilinv_c(uc_endo, d1egm, **util_para)          # Euler equation
    m_egm_agrid = (   c_egm_agrid                                               # budget constraint
                    + a_grid[np.newaxis, np.newaxis, :]) 

    max_cegm = (   m_grid[np.newaxis, np.newaxis, :]
                -  0 * d1_grid[np.newaxis, :, np.newaxis] 
                +  0 * y_grid[:, np.newaxis, np.newaxis] ) 

    max_cegm = np.maximum(max_cegm, 10 ** -6)           # impose very low non-negative consumption when upkeep too expensive         
    
    # interpolate with upper envelope, enforce borrowing limit
    Vegm, Vmegm, Vdegm, cegm, aegm = upperenv(W, Wd, m_egm_agrid, max_cegm, m_grid, a_grid, d1egm, psi, eis, eisd)
    
    # durable FOC conditional on m with durable stock in last dimension
    eval_FOC_swap = (-Vdegm + dur_price * Vmegm).swapaxes(-2,-1)
    Vegm_swap = Vegm.swapaxes(-2,-1)
    Vdegm_swap = Vdegm.swapaxes(-2,-1)
    
    coh_adj_m_swap =   (  m_grid[np.newaxis, np.newaxis, :] 
                        + dur_price * d1_grid[np.newaxis, :, np.newaxis] 
                        + np.zeros_like(Vegm)).swapaxes(-2,-1)
    dadjx, Vadjx, Vdadjx = upperenv_durable(Vegm_swap, Vdegm_swap, eval_FOC_swap, coh_adj_m_swap, d1_grid, x_grid)
    
    # interpolate onto COH grid
    coh_adj_a = coh[1,...].squeeze()

    dadj = interpolate.interpolate_y(x_grid[np.newaxis, np.newaxis, :], coh_adj_a, dadjx[:, np.newaxis, :])
    # Vadj = interpolate.interpolate_y(x_grid[np.newaxis, np.newaxis, :], coh_adj_a, Vadjx[:, np.newaxis, :])
    # Vdadj = interpolate.interpolate_y(x_grid[np.newaxis, np.newaxis, :], coh_adj_a, Vdadjx[:, np.newaxis, :])
    
    # durable choice on (adj, y, d, a) grid
    d1 = np.zeros_like(Va)

    d1[1, ...] = dadj
    d1[0, ...] = (1 - deltad_maint) * d1_grid[np.newaxis, :, np.newaxis] + np.zeros_like(dadj)
    d1 = np.maximum(d1, 10 ** - 6)
    
    d_purchase = np.zeros_like(coh)
    d_purchase[1,...] = d1[1,...]

    
    # find optimal a, c from EGM solutions given coh and d1
    n_all = np.prod(coh.shape)
    income = y_grid[np.newaxis, :, np.newaxis, np.newaxis] + np.zeros_like(coh)
    m_post_adj = coh - dur_price * d_purchase
    adjust_points_m = np.concatenate((income.reshape([n_all,1]), d1.reshape([n_all,1]), m_post_adj.reshape([n_all,1])), axis=1)
    
    a = interpn((y_grid, d1_grid, m_grid), aegm, adjust_points_m, method='linear', bounds_error=False, fill_value=None).reshape(coh.shape)
    a = np.maximum(a, a_grid[0])

    

    c = m_post_adj - a
    c = np.maximum(c, 10 ** - 6)

    # interpolation points for value functions
    all_points_a = np.concatenate((income.reshape([n_all,1]), d1.reshape([n_all,1]), a.reshape([n_all,1])), axis=1)
    # noadjust_points_a = all_points_a.reshape([2, int(n_all/2), 3])[0,...].squeeze()

    # Update value functions given c,d,a choices
    muc = margutil_c(c, d1, **util_para)
    Va = (1 + interestrate) * muc 

    Wdinterp_a = interpn((y_grid, d1_grid, a_grid), Wd, all_points_a, method='linear', bounds_error=False, fill_value=None).reshape(coh.shape)

    Vd = np.zeros_like(Va)
    Vd[1, ...] = dur_value[1, ...] * muc[1, ...]

    # Vd[1, ...] = (1 - deltad) * (1 - fc) * ( margutil_d(c[1, ...], d1[1, ...], **util_para).squeeze() 
                                # + Wdinterp_a[1, ...] ).squeeze()

    

    # Wdinterp_noadj_a = interpn((y_grid, d1_grid, a_grid), Wd, noadjust_points_a, method='linear', bounds_error=False, fill_value=None).reshape(coh.shape[1:])
    mud = margutil_d(c, d1, **util_para)
    Vd[0, ...] = ((1 - deltad_maint) * (   mud[0, ...] + Wdinterp_a[0, ...] )
                                        + dur_value[0, ...] * muc[0, ...] )
                                         
    

    Winterp_a = interpn((y_grid, d1_grid, a_grid), W, all_points_a, method='linear', bounds_error=False, fill_value=None).reshape(coh.shape)
    V = util(c, d1, **util_para) + Winterp_a

    # print(aegm.mean(), dadjx.mean() , dadj.mean(), a.mean(), (V[1,...]-V[0,...]>0).sum())
    # print((V[1,...]-V[0,...]>0).sum())
    # print((V[1,...]-V[0,...]>0).sum(axis=2).sum(axis=0))
    # print((V[1,...]-V[0,...]>0).sum(axis=1).sum(axis=0))

    # check budget constraint
    assert np.max(np.abs(a + dur_price * d_purchase + c - coh) * (c > 10 ** - 6).astype(int))<10**-8
    
    
    return V, Va, Vd, a, d1, c

def net_assets(a, d1, c, collateral_fc, borrow, deltad, d1_grid, y_grid, r, p, pbar, fc, op_cost, spread):
    # net asset position
    anet = a - collateral_fc * d1 
    aborrow = - anet * borrow

    # durable expenditure
    fixedcost = fc * np.array([0,1])[:, np.newaxis, np.newaxis, np.newaxis]
    x = d1 - (1 - deltad) * (1 - fixedcost) * d1_grid[np.newaxis, np.newaxis, :, np.newaxis]

    # budget constraint
    income = y_grid[np.newaxis, :, np.newaxis, np.newaxis]
    durexp = p * x
    nondurexp = c + pbar * op_cost * d1 + spread * aborrow 
    totexp = durexp + nondurexp
    bdgt_constr = income + r * anet - totexp

    return anet, aborrow, x, bdgt_constr, nondurexp, durexp, totexp

optad_stage = Continuous2D(backward=['Va','V','Vd'], policy=['d1','a'], f=optimal_ad,
                            name='dchoice', hetoutputs=[net_assets])   


def hh_init(coh, a_grid, d1_grid, eis, psi, eisd, collateral_fc):
    V = np.ones([2,1,1,1]) * util(0.2 * (coh[1,...] + collateral_fc * d1_grid[np.newaxis,np.newaxis,:,np.newaxis]), d1_grid[np.newaxis,np.newaxis,:,np.newaxis], psi=psi, eis=eis, eisd=eisd) / 0.01    
    Va = np.empty_like(V)
    
    Va[..., 1:-1] = (V[..., 2:] - V[..., :-2]) / (a_grid[2:] - a_grid[:-2])
    Va[..., 0] = (V[..., 1] - V[..., 0]) / (a_grid[1] - a_grid[0])
    Va[..., -1] = (V[..., -1] - V[..., -2]) / (a_grid[-1] - a_grid[-2])
    Vd = np.zeros_like(V)

    assert Va.min() > 0
    
    return V, Va, Vd

def combine_grids(min_val, max_val, split_val, n, frac):
    n_low = int(np.floor(n * frac))
    n_high = n - n_low
    grid_low = grids.agrid(split_val, n_low, min_val)
    grid_high = np.linspace(split_val, max_val, n_high + 1)[1:]
    grid = np.concatenate((grid_low,grid_high), axis=0)
    return grid

def make_grids(rho_z, sd_z, n_z, min_a, max_a, n_a, min_d, max_d, n_d, min_m, max_m, n_m, n_x):
    z_grid, z_dist, z_markov = grids.markov_rouwenhorst(rho_z, sd_z, n_z)
    a_grid = grids.agrid(max_a, n_a, min_a)
    # a_grid_low = np.array([min_a])
    # a_grid_high = grids.agrid(max_a, n_a - 1, 10 ** -4)
    # a_grid = np.concatenate((a_grid_low,a_grid_high), axis=0)

    m_grid = grids.agrid(max_m, n_m, min_m)
    d1_grid = grids.agrid(max_d, n_d, min_d)
    x_grid = grids.agrid(max_d + max_a, n_x, min_a + min_d)

    adjust_grid = np.linspace(0, 1, 2)
    return z_grid, z_dist, z_markov, a_grid, m_grid, d1_grid, x_grid, adjust_grid


def labor_income(a_grid, m_grid, d1_grid, z_grid, r, W, p, deltad, pbar, fc, maint_cost, op_cost, collateral, spread, Tr):
    y_grid = z_grid * W + Tr

    y = y_grid[np.newaxis, :, np.newaxis, np.newaxis]           # on (adj, z, d, a)

    collateral_fc = collateral * (1 - deltad) * (1 - fc) * pbar

    # effective interest rate based on whether household was borrowing
    atilde_on_grid = a_grid[np.newaxis, np.newaxis, np.newaxis, :] - collateral_fc * d1_grid[np.newaxis, np.newaxis, :, np.newaxis]
    borrow = (atilde_on_grid < 0).astype(int)
    if  spread>0:
        interestrate = r  + spread * (atilde_on_grid < 0)
    else:
        interestrate = r
    
    # durable values
    deltad_maint = deltad * (1 - maint_cost)
    maint_cost_ratio = p * maint_cost * deltad / (1 - deltad_maint)

    dur_maint = pbar * op_cost + maint_cost_ratio
    dur_price = p + pbar * op_cost  - collateral_fc

    # cash on hand after durable upkeep and durable sales
    dur_value = np.zeros([2,1,len(d1_grid),len(a_grid)])
    dur_value[1, ...] = p * (1 - deltad) * (1 - fc)  - (1 + interestrate) * collateral_fc
    dur_value[0, ...] = - dur_maint * (1 - deltad_maint) - collateral_fc * (interestrate + deltad_maint)

    coh = (1 + interestrate) * a_grid[np.newaxis, np.newaxis, np.newaxis, :] + dur_value * d1_grid[np.newaxis, np.newaxis, :, np.newaxis] + y
    
    return y, y_grid, coh, dur_value, deltad_maint, dur_price, interestrate, borrow, collateral_fc

@simple
def market_clearing(B_ss, ANET, ABORROW, r, W, Tr, C, p, D1, X, op_cost, pbar, spread, collateral): #, D, p
    asset_mkt = ANET - B_ss 
    borrow_cost = spread * ABORROW
    bdgt_constr2 = W + Tr + r * ANET - C - p * X - pbar * op_cost * D1 - borrow_cost
    
    return asset_mkt, bdgt_constr2, borrow_cost

hh = StageBlock([prod_stage, adjust_stage, optad_stage], name='hh',
                backward_init=hh_init, hetinputs=[make_grids, labor_income])

cali = dict()

rhoYq = 0.966
rhoYm = 0.966**(1/3)
sigmaYq = 0.5
sigmaYm = sigmaYq * 3 / np.sqrt(1 + (1 + rhoYq) ** 2 + (1 + rhoYq + rhoYq**2) ** 2 + rhoYq ** 2 * (1 + rhoYq) ** 2  + rhoYq ** 2)

cali = {'taste_shock': 1E-3 * 1000 / 5, 'r': 0.02/12,
        'eis': 1, 'eisd': 1, 
        'psi': 0.11, 'deltad': 0.015, 'fc': 0.045,
        'collateral': 0.8, 'spread': 0.02 / 12,
        'maint_cost': 0.466, 'op_cost': 0.018, 
        'Tr': 0,
        # 'beta': 0.9948110014496739,
        'beta': 0.995,
        # 'B_ss': 0.1,
        'rho_z': rhoYm, 'sd_z': sigmaYm, 'n_z': 2, 
        'min_a': 0.0, 'max_a': 120/12, 'n_a': 150,
        'min_m': 0.1/12, 'max_m': 120/12, 'n_m': 400,
        'min_d': 0.01, 'max_d': 12/12, 'n_d': 50,
        'n_x': 400,
        'W': 1/12, 'p': 1, 'pbar': 1}             

# fc_model = create_model([hh, market_clearing])
# unknowns_ss = {'beta': (0.987/(1+cali['r']), 0.995/(1+cali['r']))}
# targets_ss = {'asset_mkt': 0.}
# ss = fc_model.solve_steady_state(cali, unknowns_ss, targets_ss, solver='brentq')
# assert np.abs(ss['bdgt_constr2'])<10 ** - 5
# jac = fc_model.solve_jacobian(ss, ['Y'], ['bdgt_constr2'], ['r'], T=5)


ss = hh.steady_state(cali)
print(f"Aggregate durables: {ss['D1']/ss['NONDUREXP']/12:0.2f}")
print(f"Aggregate assets: {ss['ANET']:0.2f}")
print(f"Aggregate BC: {ss['BDGT_CONSTR']:0.4f}")

jac = hh.jacobian(ss, inputs=['Tr','p','r'], T=30)

print(jac['TOTEXP']['Tr'][0:3,0].sum(), jac['NONDUREXP']['Tr'][0:3,0].sum(), jac['DUREXP']['Tr'][0:3,0].sum())
print(jac['DUREXP']['p'][0:6,6:30].sum() / (ss['DUREXP'] * 6) )

def fig1(ss, dmax=40, dmin=0, yz=0, az=0, figsize=0.6):
    a_grid = ss.internals['hh']['a_grid']
    m_grid = ss.internals['hh']['m_grid']
    d1_grid = ss.internals['hh']['d1_grid']

    a = ss.internals['hh']['dchoice']['a']
    da = a - a_grid
    c = ss.internals['hh']['dchoice']['c']
    d1 = ss.internals['hh']['dchoice']['d1']
    V = ss.internals['hh']['dchoice']['V']
    Vd = ss.internals['hh']['dchoice']['Vd']
    Va = ss.internals['hh']['dchoice']['Va']
    P = ss.internals['hh']['adjust']['law_of_motion'].P
        # FOC = ss.internals['hh']['dchoice']['eval_FOC_swap']
        

    fig, axes = plt.subplots(2, 2, figsize=(8*figsize, 8*figsize))
    ax = axes.flatten()
    
    label = [r'fixed cost model']
    adjlabels = [r'adjusting', r'not adjusting'] 

    for adj, l in zip([1, 0], adjlabels):
        ax[0].plot(d1_grid[dmin:dmax], a[adj, yz, dmin:dmax, az], label=l, linewidth=2)
    # ax[0].plot(d1_grid[dmin:dmax], d1_grid[dmin:dmax], color='gray', linestyle=':')
    ax[0].legend(frameon=False)
    ax[0].axhline(0, color='gray', linestyle=':')
    ax[0].set_title('Assets')

    # for i, ltemp in zip(models, labels):
    #     for adj, l in zip([1, 0], adjlabels):
    #         ax[1].plot(a_grid[:amax], c[adj, yz, dz, :amax], label=l, linewidth=2)
    # ax[1].set_title('Consumption')
    ax[1].plot(d1_grid[dmin:dmax], Vd[1, yz, dmin:dmax, az], label='Vd', linewidth=2)
    ax[1].plot(d1_grid[dmin:dmax], Va[1, yz, dmin:dmax, az], label='Va', linewidth=2)
    # ax[1].plot(d1_grid[dmin:dmax], V[1, yz, dmin:dmax, az], label='V', linewidth=2)
    ax[1].set_title('Durable Value')
    ax[1].legend(frameon=False)

    for adj, l in zip([1, 0], adjlabels):
        ax[2].plot(d1_grid[dmin:dmax], d1[adj, yz, dmin:dmax, az], label=l, linewidth=2)
    ax[2].set_title('Durable')
    ax[2].legend(frameon=False)

    # for i, l in zip(models, labels):
    #     user_cost = ss[i]['p'] * (ss[i]['r'] + ss[i]['deltad']) / (1 + ss[i]['r']) #pbar * op_cost
    #     opt_ratio = user_cost / ss[i]['psi']
    #     ax[3].plot(a_grid[:amax], c[i][1, yz, dz, :amax] / d1[i][1, yz, dz, :amax], label=l, linewidth=2)
    #     ax[3].plot(a_grid[:amax], opt_ratio * np.ones_like(c[i][1, yz, dz, :amax]), label='frictionless', linewidth=2)
    # ax[3].set_title('Nondurable / Durable')
    # ax[3].legend(frameon=False)

    ax[3].plot(d1_grid[dmin:dmax], P[1, 0, yz, dmin:dmax, az], label=l, linewidth=2)
    ax[3].set_title('Adjustment Probability')
    
    for k in ax:
        k.set_xlabel('Durables')

    plt.tight_layout()
    plt.show()

    print((P[1,...].squeeze() * ss.internals['hh']['dchoice']['D']).sum())


#     # fig, ax = plt.subplots()
#     # ax.plot(d1_grid[dmin:dmax], FOC[i][yz, az, dmin:dmax], label=l, linewidth=2)
#     # ax.set_title('FOC')
#     # ax.set_xlabel('Durables')



fig1(ss, dmax=cali['n_d']-1, yz=0, az=0, figsize=0.8) 
# fig1(ss, dmax=cali['n_d']-1, yz=0, az=int(cali['n_a']/2), figsize=0.8)    
# fig1(ss, dmax=cali['n_d']-1, yz=cali['n_z'] - 1, az=cali['n_a'] - 1, figsize=0.8)    
# # %%

# %%
