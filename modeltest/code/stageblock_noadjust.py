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
def util(c, d, psi=0, eis=1):
    u = np.log(c) + psi * np.log(d)
    return u

@njit
def margutil_c(c, d, psi=0, eis=1):
    muc = c ** -1
    return muc  

@njit
def margutilinv_c(muc, d, psi=0, eis=1):
    c = muc ** -1
    return c       

@njit
def margutil_d(c, d, psi=0, eis=1):
    mud = psi * d ** -1
    return mud        

@njit
def upperenv_vec(W, Wd, m_endo, max_c, m_grid, a_grid, d, *args):
    
    """Interpolate value function and consumption to exogenous grid."""
    n_b, n_m = max_c.shape
    n_b2, n_a = W.shape

    assert max_c.min()>0

    assert np.abs(n_b-n_b2)<10**-8

    a = np.zeros_like(max_c)
    c = np.zeros_like(max_c)
    V = -np.inf * np.ones_like(max_c)
    Vd = -np.inf * np.ones_like(max_c)
    Vm = -np.inf * np.ones_like(max_c)

    # loop over other states, collapsed into single axis
    for ib in range(n_b):
        d0 = d[ib, 0]

        # loop over segments of endogenous asset grid from EGM (not necessarily increasing)
        for ja in range(n_a - 1):
            m_low, m_high = m_endo[ib, ja], m_endo[ib, ja + 1]
            W_low, W_high = W[ib, ja], W[ib, ja + 1]
            Wd_low, Wd_high = Wd[ib, ja], Wd[ib, ja + 1]
            ap_low, ap_high = a_grid[ja], a_grid[ja + 1]
            
           # loop over exogenous asset grid (increasing) 
            for im in range(n_m):  
                mcur = m_grid[im]
                max_c_cur = max_c[ib, im]
                
                interp = (m_low <= mcur <= m_high) 
                extrap = (ja == n_a - 2) and (mcur > m_endo[ib, n_a - 1])

                # exploit that a_grid is increasing
                if (m_high < mcur < m_endo[ib, n_m - 1]):
                    break

                if interp or extrap:
                    W0 = interpolate.interpolate_point(mcur, m_low, m_high, W_low, W_high)
                    Wd0 = interpolate.interpolate_point(mcur, m_low, m_high, Wd_low, Wd_high)
                    a0 = interpolate.interpolate_point(mcur, m_low, m_high, ap_low, ap_high)
                    c0 = max_c_cur - a0
                    V0 = util(c0, d0, *args) + W0
                    Vm0 = margutil_c(c0, d0, *args)
                    Vd0 = margutil_d(c0, d0, *args) + Wd0

                    # upper envelope, update if new is better
                    if V0 > V[ib, im]:
                        a[ib, im] = a0 
                        c[ib, im] = c0
                        V[ib, im] = V0
                        Vm[ib, im] = Vm0
                        Vd[ib, im] = Vd0

        # Enforce borrowing constraint
        im = 0
        while im < n_m and m_grid[im] <= m_endo[ib, 0]:
            a[ib, im] = a_grid[0]
            c[ib, im] = max_c[ib, im]
            V[ib, im] = util(c[ib, im], d0, *args) + W[ib, 0]
            Vm[ib, im] = margutil_c(c[ib, im], d0, *args)
            Vd[ib, im] = margutil_d(c[ib, im], d0, *args) + Wd[ib, 0]
            im += 1

    assert c.min()>0
    assert a.min()>a_grid[0] - 10 ** - 10

    return V, Vm, Vd, c, a                                               


def upperenv(W, Wd, m_endo, max_c, m_grid, a_grid, d, *args):
    # collapse (z, d, a) into (b, a)
    shapea = W.shape
    W = W.reshape((-1, shapea[-1]))
    Wd = Wd.reshape((-1, shapea[-1]))
    m_endo = m_endo.reshape((-1, shapea[-1]))
    d = d.reshape((-1, shapea[-1]))

    shapem = max_c.shape
    max_c = max_c.reshape((-1, shapem[-1]))
    
    
    V, Vm, Vd, c, a = upperenv_vec(W, Wd, m_endo, max_c, m_grid, a_grid, d, *args)
    
    # report on (z, d, a)
    return V.reshape(shapem), Vm.reshape(shapem), Vd.reshape(shapem), c.reshape(shapem), a.reshape(shapem)                   


@njit
def upperenv_durable_vec(V_FOC, eval_FOC, coh, d1_grid, x_grid, d):

    """Interpolate durable expenditure to exogenous grid."""
    n_b, n_m, n_d = eval_FOC.shape
    n_x = x_grid.shape[0]
    n_d = d1_grid.shape[0]


    # assert V_FOC.shape == eval_FOC.shape
    # assert coh.shape == eval_FOC.shape
    assert coh.min()>0

    
    # d = np.zeros([n_b,n_x])
    # V = -np.inf * np.ones([n_b,n_x])
    V = -np.inf * np.ones_like(d)
    
    # loop over other states, collapsed into single axis
    for ib in range(n_b):

        xlist = []
        dlist = []
        Vlist = []

        for id in range(n_m):
            # m0 = m_grid[id]

            # loop over segments of FOC (not necessarily increasing)
            for jd in range(n_d - 1):
                foc_low, foc_high = eval_FOC[ib, id, jd], eval_FOC[ib, id, jd + 1]

                if np.sign(foc_low) != np.sign(foc_high):
                    x_low, x_high = coh[ib, id, jd], coh[ib, id, jd + 1]
                    d_low, d_high = d1_grid[jd], d1_grid[jd + 1]
                    V_FOC_low, V_FOC_high = V_FOC[ib, id, jd], V_FOC[ib, id, jd + 1]

                    x0 = interpolate.interpolate_point(0, foc_low, foc_high, x_low, x_high)
                    d0 = interpolate.interpolate_point(0, foc_low, foc_high, d_low, d_high)
                    V0 = interpolate.interpolate_point(0, foc_low, foc_high, V_FOC_low, V_FOC_high)

                    dlist.append(d0)
                    xlist.append(x0)
                    Vlist.append(V0)

        # now have combos given coh0, choice (m0, d0) satisfies FOC
        # with unique solution, will just have one d0 for each x0 and increasing in x0
        n_jm = len(xlist)
        
        # loop over adjecent (x,d) solutions (not necessarily increasing)
        for jm in range(n_jm - 1):
            x_low, x_high = xlist[jm], xlist[jm + 1]
            d_low, d_high = dlist[jm], dlist[jm + 1]
            V_low, V_high = Vlist[jm], Vlist[jm + 1]

            # loop over exogenous cash on hand grid (increasing) 
            for ix in range(n_x):  
                xcur = x_grid[ix]
                
                interp = (x_low <= xcur <= x_high)
                extrap = (jm == n_jm - 2) and (xcur > xlist[n_jm - 1])

                # exploit that x_grid is increasing
                if (x_high < xcur < xlist[n_jm - 1]):
                    break

                if interp or extrap:
                    d1 = interpolate.interpolate_point(xcur, x_low, x_high, d_low, d_high)
                    V1 = interpolate.interpolate_point(xcur, x_low, x_high, V_low, V_high)

                    # upper envelope, update if new is better
                    if V1 > V[ib, ix]:
                        d[ib, ix] = d1
                        V[ib, ix] = V1

            # Extrapolate at boundary
            ix = 0
            xmin = xlist[0]
            while ix < n_x and x_grid[ix] <= xmin:
                d[ib, ix] = d1_grid[0]
                ix += 1 

            if ix >= 1:
                d[ib, ix-1] = interpolate.interpolate_point(x_grid[ix-1], xlist[0], xlist[1], dlist[0], dlist[1])      
    
    return d, V

def upperenv_durable(V_FOC, eval_FOC, coh, d1_grid, x_grid):
    # collapse (z, m, d) into (b, m, d)
    shapem = V_FOC.shape
    V_FOC = V_FOC.reshape((-1, shapem[-2], shapem[-1]))
    eval_FOC = eval_FOC.reshape((-1, shapem[-2], shapem[-1]))
    coh = coh.reshape((-1, shapem[-2], shapem[-1]))

    shapex = (np.prod(shapem[:-2]), x_grid.shape[0])
    shapeout = shapem[:-2] + (x_grid.shape[0],)

    assert shapem[-1] == d1_grid.shape[0]
    assert V_FOC.shape[0] == np.prod(shapem[:-2])

    d = np.zeros(shapex)
    
    d, V = upperenv_durable_vec(V_FOC, eval_FOC, coh, d1_grid, x_grid, d)
    
    # report on (z, x) 
    return d.reshape(shapeout), V.reshape(shapeout)   


def optimal_ad(V, Va, Vd, coh, interestrate, adj_dur_value, dur_maint, dur_price, beta, psi, eis, a_grid, m_grid, d1_grid, y_grid, x_grid):

    util_para = {'psi': psi, 'eis': eis}

    # Continuation value function
    W = beta * V                                                  # end-of-stage vfun
    Wd = beta * Vd                                                 # end-of-stage vfun
    uc_endo = beta * Va                                           # envelope condition

    d1 = d1_grid[np.newaxis,:,np.newaxis] + np.zeros_like(V)
    
    c_endo = margutilinv_c(uc_endo, d1, **util_para)                # Euler equation
    m_endo = ( c_endo                                               # budget constraint
             + a_grid[np.newaxis, np.newaxis, :] 
             + dur_maint * d1_grid[np.newaxis, :, np.newaxis]) 
    
    max_cegm = (   m_grid[np.newaxis, np.newaxis, :]
                -  dur_maint * d1_grid[np.newaxis, :, np.newaxis]
                +  0 * y_grid[:, np.newaxis, np.newaxis] )

    max_cegm = np.maximum(max_cegm, 10 ** -6)           # impose very low non-negative consumption when upkeep too expensive                          

    # interpolate with upper envelope, enforce borrowing limit
    Vegm, Vmegm, Vdegm, cegm, aegm = upperenv(W, Wd, m_endo, max_cegm, m_grid, a_grid, d1, psi, eis)
    

    # test code
    uc_endo = beta * Va 
    c_endo = margutilinv_c(uc_endo, d1, **util_para)
    exp_endo =  ( c_endo 
                + a_grid[np.newaxis,np.newaxis,:] 
                + dur_maint * d1_grid[np.newaxis, :, np.newaxis]) 
    aegm2 = interpolate.interpolate_y(exp_endo, m_grid[np.newaxis, np.newaxis, :], a_grid)    
    aegm2 = np.maximum(aegm2, a_grid[0])
    assert np.max(np.abs(aegm - aegm2))<10 ** -8

    # durable FOC conditional on m with durable stock in last dimension
    # Vdegm = u_d + beta W_d = (p + \nu - \lambda * \bar{p}) * u_c = (p + \nu - \lambda * \bar{p}) V_m
    eval_FOC = (-Vdegm + dur_price * Vmegm).swapaxes(-2,-1)
    degm = interpolate.interpolate_y(eval_FOC, np.zeros(eval_FOC.shape[:-1]+(1,)), d1_grid[np.newaxis, np.newaxis, :])   
    degm = degm.squeeze()


    # interpolate onto coh grid
    coh_endo = m_grid[np.newaxis, :] + (dur_price - dur_maint) * degm
    d1 = interpolate.interpolate_y(coh_endo[:, np.newaxis, :], coh, degm[:, np.newaxis, :])
    d1interp = np.maximum(d1, d1_grid[0])

    d2 = interpolate.interpolate_y(coh_endo, x_grid[np.newaxis, :], degm)
    d2 = np.maximum(d2, d1_grid[0])

    coh_adj =   ( m_grid[np.newaxis, np.newaxis, :] 
                + (dur_price - dur_maint) * d1_grid[np.newaxis, :, np.newaxis] 
                + np.zeros_like(Vegm)).swapaxes(-2,-1)
    dadj, Vadj = upperenv_durable(Vegm.swapaxes(-2,-1), eval_FOC, coh_adj, d1_grid, x_grid)
    
    d1 = interpolate.interpolate_y(x_grid[np.newaxis, np.newaxis, :], coh, dadj[:, np.newaxis, :])
    
    # look for small difference between interpolation and upper envelope method
    # print(np.max(np.abs(d1interp - d1)))
    # print(np.max(np.abs(dstar - d2)))
    # print(np.max(np.abs(dstar - dadj)))
    # print(dadj[0,:10],dstar[0,:10])
    # print(weight[0,:10], dadj_loc_c[0,:10])
    # print(dadj[0,-5:],dstar[0,-5:])
    # assert np.max(np.abs(d2 - dadj))< 10 **-2
    # assert np.max(np.abs(d1interp - d1))< 10 **-2
    
    # find optimal a, c from EGM solutions
    n_all = np.prod(coh.shape)
    income = y_grid[:, np.newaxis, np.newaxis] + np.zeros_like(coh)
    m_endo = coh - (dur_price - dur_maint) * d1
    adjust_points = np.concatenate((income.reshape([n_all,1]), d1.reshape([n_all,1]), m_endo.reshape([n_all,1])), axis=1)
    a = interpn((y_grid, d1_grid, m_grid), aegm, adjust_points, method='linear', bounds_error=False, fill_value=None).reshape(coh.shape)
    a = np.maximum(a, a_grid[0])
    c = m_endo - a - dur_maint * d1
    c = np.maximum(c, 10 ** - 6)


    Va = (1 + interestrate) * margutil_c(c, d1, psi, eis)  
    Vd = adj_dur_value * margutil_c(c, d1, psi, eis) 


    adjust_points = np.concatenate((income.reshape([n_all,1]), d1.reshape([n_all,1]), a.reshape([n_all,1])), axis=1)
    Vprime = interpn((y_grid, d1_grid, a_grid), W, adjust_points, method='linear', bounds_error=False, fill_value=None).reshape(coh.shape)
    
    V = util(c, d1, **util_para) + Vprime

    # check budget constraint
    assert np.max(np.abs(a + dur_price * d1 + c - coh))<10**-8
    
    return V, Va, Vd, a, d1, c

def net_assets(a, d1, collateral, pbar, borrow):
    # net asset position
    anet = a - collateral * pbar * d1 
    aborrow = anet * borrow
    return anet, aborrow

optad_stage = Continuous2D(backward=['Va','V','Vd'], policy=['d1','a'], f=optimal_ad,
                            name='dchoice', hetoutputs=[net_assets])                              


def hh_init(coh, a_grid, d1_grid, eis, psi):
    V = util(0.2 * coh, d1_grid[np.newaxis,:,np.newaxis], eis, psi) / 0.01    
    Va = np.empty_like(V)
    
    Va[..., 1:-1] = (V[..., 2:] - V[..., :-2]) / (a_grid[2:] - a_grid[:-2])
    Va[..., 0] = (V[..., 1] - V[..., 0]) / (a_grid[1] - a_grid[0])
    Va[..., -1] = (V[..., -1] - V[..., -2]) / (a_grid[-1] - a_grid[-2])
    Vd = np.zeros_like(V)
    
    return V, Va, Vd

def make_grids(rho_z, sd_z, n_z, min_a, max_a, n_a, min_d, max_d, n_d, min_m, max_m, n_m, n_x):
    z_grid, z_dist, z_markov = grids.markov_rouwenhorst(rho_z, sd_z, n_z)
    a_grid = grids.agrid(max_a, n_a, min_a)
    m_grid = grids.agrid(max_m, n_m, min_m)
    d1_grid = grids.agrid(max_d, n_d, min_d)
    x_grid = grids.agrid(max_d + max_a, n_x, min_a + min_d)
    return z_grid, z_dist, z_markov, a_grid, m_grid, d1_grid, x_grid


def labor_income(a_grid, m_grid, d1_grid, z_grid, deltad, r, w, p, pbar, maint_cost, op_cost, collateral, spread):
    y_grid = z_grid * w 

    y = y_grid[:, np.newaxis, np.newaxis]           # on (adj, z, d, a)
    coh_m_only = m_grid[ np.newaxis, np.newaxis, :] + 0 * d1_grid[ np.newaxis, :, np.newaxis] + 0 * y

    # effective interest rate based on whether household was borrowing
    atilde_on_grid = a_grid[np.newaxis, np.newaxis, :] - collateral * pbar * d1_grid[ np.newaxis, :, np.newaxis]
    borrow = (atilde_on_grid < 0).astype(int)
    if  spread>0:
        interestrate = r  + spread * (atilde_on_grid < 0)
    else:
        interestrate = r

    deltad_maint = deltad * (1 - maint_cost)
    maint_cost_ratio = p * maint_cost * deltad / (1 - deltad_maint)

    dur_maint = pbar * op_cost + maint_cost_ratio
    dur_price = p + pbar * op_cost  - collateral * pbar

    noadj_dur_value = dur_maint + collateral * pbar * (interestrate + deltad_maint) / ((1 - deltad_maint))
    adj_dur_value = p * (1 - deltad)  - (1 + interestrate) * collateral * pbar

    # assert dur_upkeep * d1_grid[-1] < m_grid[0], 'not enough cash on hand'

    coh = (1 + interestrate) * a_grid[ np.newaxis, np.newaxis, :] + adj_dur_value * d1_grid[ np.newaxis, :, np.newaxis] + y 

    return y, y_grid, coh_m_only, coh, adj_dur_value, noadj_dur_value, deltad_maint, dur_maint, dur_price, interestrate, borrow



hh = StageBlock([prod_stage, optad_stage], name='hh',
                backward_init=hh_init, hetinputs=[make_grids, labor_income])

print(hh)
print(f"Inputs: {hh.inputs}")
print(f"Outputs: {hh.outputs}")

cali = dict()

cali['ck'] = {'taste_shock': 1E-3, 'r': 0.02/12, 'beta': 0.975, 'eis': 1, 
               'psi': 0.1, 'deltad': 0.015, 'fc': 0,
               'collateral': 0.8, 'spread': 0,
               'maint_cost': 0.466, 'op_cost': 0.018, 
               'rho_z': 0.95, 'sd_z': 0.5, 'n_z': 2, 
               'min_a': 0.0, 'max_a': 40, 'n_a': 50,
               'min_m': 0.1, 'max_m': 40, 'n_m': 800,
               'min_d': 0.01, 'max_d': 24, 'n_d': 50,
               'n_x': 800,
               'w': 1.0, 'p': 1, 'pbar': 1}
# cali['ck']['min_m'] = cali['ck']['max_d'] * (cali['ck']['op_cost'] + cali['ck']['op_cost'])

# cali['ck'] = deepcopy(cali['sim'])
# cali['ck']['vphi'] = 1.0

# cali['ck_smooth'] = deepcopy(cali['ck'])
# cali['ck_smooth']['taste_shock'] = 1E-2

ss, jac = dict(), dict()
for i in cali.keys():
    ss[i] = hh.steady_state(cali[i])
    print(f"Aggregate durables in model {i}: {ss[i]['D1']:0.2f}")
    print(f"Aggregate assets in model {i}: {ss[i]['ANET']:0.2f}")
    jac['ck'] = hh.jacobian(ss['ck'], inputs=['r'], T=5)

def fig1(ss, amax=150, amin=0, yz=0, dz=0, figsize=0.6):
    a_grid = ss['ck'].internals['hh']['a_grid']
    a, da, c, d1, V = dict(), dict(), dict(), dict(), dict()
    models = ['ck']
    for i in models:
        a[i] = ss[i].internals['hh']['dchoice']['a']
        da[i] = a[i] - a_grid
        c[i] = ss[i].internals['hh']['dchoice']['c']
        d1[i] = ss[i].internals['hh']['dchoice']['d1']
        V[i] = ss[i].internals['hh']['dchoice']['V']

    fig, axes = plt.subplots(2, 2, figsize=(8*figsize, 8*figsize))
    ax = axes.flatten()
    
    labels = [r'no adjustment cost'] 

    for i, l in zip(models, labels):
        ax[0].plot(a_grid[:amax], a[i][yz, dz, :amax], label=l, linewidth=2)
    ax[0].plot(a_grid[:amax], a_grid[:amax], color='gray', linestyle=':')
    ax[0].legend(frameon=False)
    ax[0].axhline(0, color='gray', linestyle=':')
    ax[0].set_title('Assets')

    for i, l in zip(models, labels):
        ax[1].plot(a_grid[:amax], c[i][yz, dz, :amax], label=l, linewidth=2)
    ax[1].set_title('Consumption')

    for i, l in zip(models, labels):
        ax[2].plot(a_grid[:amax], d1[i][yz, dz, :amax], label=l, linewidth=2)
    ax[2].set_title('Durable')

    for i, l in zip(models, labels):
        user_cost = ss[i]['p'] * (ss[i]['r'] + ss[i]['deltad']) / (1 + ss[i]['r']) + ss[i]['pbar'] * ss[i]['op_cost']
        opt_ratio = user_cost / ss[i]['psi']
        ax[3].plot(a_grid[:amax], c[i][yz, dz, :amax] / d1[i][yz, dz, :amax], label=l, linewidth=2)
        ax[3].plot(a_grid[:amax], opt_ratio * np.ones_like(c[i][yz, dz, :amax]), label='frictionless', linewidth=2)
    ax[3].set_title('Nondurable / Durable')
    ax[3].legend(frameon=False)
    
    for k in ax:
        k.set_xlabel('assets')

    plt.tight_layout()
    plt.show()

fig1(ss, amax=85, yz=1, dz=0, figsize=0.8)   
# %%
