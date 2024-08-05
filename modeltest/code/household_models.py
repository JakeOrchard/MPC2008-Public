import numpy as np

import sequence_jacobian.utilities as utils
from newton_solver2 import newton_solver2
from sequence_jacobian import simple, het, create_model, estimation, interpolate, misc
from sequence_jacobian.utilities.interpolate import interpolate_y

from scipy.interpolate import interpn
from customutils import setmin3

from utilityfunction import f_u, f_invuc, f_uc, f_ud, f_invuc_dc, f_uc_dc

def household_init(a_grid, d1_grid, x_grid, y, r, p, eis, psi, deltad):
    coh = (1 + r) * a_grid[np.newaxis, np.newaxis, :] + p * (1 - deltad) * d1_grid[np.newaxis, :, np.newaxis] + y[:, np.newaxis, np.newaxis]
    Va = (1 + r) * (0.1 * psi * coh) ** (-1 / eis)
    Vd = (1-deltad) * Va
    Vx = Va / (1+r)
    
    V = np.cumsum(Va, axis=2) + np.cumsum(Vd, axis=1)
    
    return Vx, Va, Vd, V

@het(exogenous='Pi', policy=['d1','a'], backward=['V','Vd','Va','Vx'], backward_init=household_init)
def household_nofc(Vx_p, Va_p, Vd_p, V_p,  a_grid, d1_grid, s_grid, x_grid, y, r, r_p, spread, op_cost, collateral, pbar, p, p_p, user_cost, user_cost_borrow, beta, eis, psi, deltad, xi, d1_mesh, a_mesh, y_mesh, y_mesh_yx, x_mesh_yx, nE, nA, nX, nD, dmin):

    nN = np.prod(y_mesh.shape)
    nNyx = np.prod(y_mesh_yx.shape)

    # effective interest rate based on whether household was borrowing
    atilde_mesh = a_mesh - collateral * pbar * d1_mesh

    if  spread>0:
        borrowed = (atilde_mesh<0)
        interestrate0 = r * np.ones(atilde_mesh.shape)
        interestrate0[borrowed] = r + spread
    else:
        interestrate0 = r

    # util_params = {'eis': eis, 'xi': xi, 'psi': psi}
    effective_dur_price = p + pbar * op_cost - collateral * pbar 
    # effective_dur_price_borrow = p + pbar * op_cost - collateral * pbar / (1 + r + spread)
    effective_dur_price_kink = p + pbar * op_cost
    
    # future_dur_value = (p_p * (1 - deltad) - collateral * pbar) 
    future_dur_value = (p_p * (1 - deltad) - (1 + r_p) * collateral * pbar) 
    future_dur_value_borrow = (p_p * (1 - deltad) - (1 + r_p + spread) * collateral * pbar) 

    current_dur_value = (p * (1 - deltad) - (1 + interestrate0) * collateral * pbar)

    x_adjust_points = np.concatenate((y_mesh_yx.reshape([nNyx,1]), 
                                      np.zeros([nNyx,1]), 
                                      x_mesh_yx.reshape([nNyx,1])), 
                                      axis=1)

    Vx_p_interp = interpn((y, d1_grid, a_grid), Vx_p, x_adjust_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nX])
    Vdx_p = interpn((y, d1_grid, a_grid), Vd_p, x_adjust_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nX])
    
    # === STEP 1: Find unconstraind solution ===
    # === (a) Find EGM solution assuming no borrowing===
    # (take discounted expectation of tomorrow's value function)
    # c' given d', a'
    def solve_egm_x(Vx_p = Vx_p_interp, 
                    beta = beta, 
                    interestrate_p = r_p, 
                    user_cost = user_cost,
                    dur_price = effective_dur_price, 
                    dur_value_p = future_dur_value,
                    savegrid = x_grid[np.newaxis, :], 
                    cohgrid = x_grid[np.newaxis, :],
                    ):

        c_prime = psi / (beta * (1+interestrate_p) * Vx_p)
        d1_prime = (1-psi) / psi * c_prime / user_cost

        # x' = (1 + r') (a' - lambda * pbar * d') + (1 - delta) p'd' 
        a_prime = savegrid - dur_value_p * d1_prime

        coh0 = (a_prime - y[:, np.newaxis]
            + c_prime + dur_price * d1_prime )

        aegm = interpolate_y(coh0, cohgrid, a_prime)
        aegm = np.maximum(aegm, a_grid.min())

        return aegm

    astar = solve_egm_x()
    
    # === STEP 1(b): Find solution for c',d' given a' ===
    # how to divide c',d' when a' is given
    def solve_foc_x(astar = astar,
                    Vx_p = Vx_p_interp,
                    coh_noy = x_grid[np.newaxis,:],
                    income = y[:,np.newaxis],
                    candidate_grid = s_grid[np.newaxis,np.newaxis,:],
                    beta = beta,
                    interestrate_p = r_p,
                    dur_price = effective_dur_price,
                    dur_value_p = future_dur_value,
                    ):
                     

        spending = (coh_noy + income - astar)
    
        dcandidates = spending[:,:,np.newaxis] / dur_price * candidate_grid
        ccandidates = spending[:,:,np.newaxis] * (1 - candidate_grid)

        Vx_p_bc = interpolate_y(coh_noy[np.newaxis,:,:], 
                                (1 + interestrate_p) * astar[:,:,np.newaxis] + dur_value_p * dcandidates, 
                                Vx_p[:,np.newaxis,:])                            

        evalfoc = -(1 - psi) / dcandidates + dur_price * psi / ccandidates - beta * dur_value_p * Vx_p_bc
        # evalfoc = -(1 - psi) / dcandidates + dur_price * psi/ccandidates - beta * dur_value_p * Vx_p_bc
        # evalfoc = -(1 - psi) / dcandidates + dur_price * psi/ccandidates - dur_value_p / (1+r) * psi/ccandidates
        
        d_bc = np.squeeze(utils.interpolate.interpolate_y(evalfoc, np.zeros([nE,nX,1]), dcandidates))
        c_bc = spending - dur_price * d_bc

        return d_bc, c_bc

    d_bc, c_bc = solve_foc_x()
    
    # atilde_bc = astar - collateral * pbar * d_bc

    # === STEP 1(c): Interpolate onto a,d grid ===
    # interpolate onto a,d grid:    
    def interp_to_base_grid(d_coh = d_bc,
                            c_coh = c_bc,
                            coh_noy = (1 + interestrate0) * a_mesh + current_dur_value * d1_mesh,
                            income = y_mesh,
                            dur_price = effective_dur_price,
                            borrowing_d = collateral * pbar
                            ):
        x_adjust_points = np.concatenate((income.reshape([nN,1]), coh_noy.reshape([nN,1])), axis=1)

        d1 = interpn((y, x_grid), d_coh, x_adjust_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nD,nA])
        c = interpn((y, x_grid), c_coh, x_adjust_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nD,nA])

        a = coh_noy + income - c - dur_price * d1
        atilde = a - borrowing_d * d1

        return c, d1, a, atilde

    c, d1, a, atilde = interp_to_base_grid()
    assert a.min()>-10**-8, 'error ' + np.array2string(a.min())
    
    if spread>0:
        # === STEP 2: Solution assuming borrowing with spread > 0 ===
        # === (a) EGM part of algorithm ===
        astar_borrow = solve_egm_x(interestrate_p = r_p + spread,  
                                   user_cost = user_cost_borrow,
                                   dur_value_p = future_dur_value_borrow)

        # === STEP 2(b): Find solution for c',d' given a' assuming borrowing ===
        # how to divide c',d' when a' is given
        d_bc_borrow, c_bc_borrow = solve_foc_x(astar = astar_borrow,
                                               interestrate_p = r_p + spread,
                                               dur_value_p = future_dur_value_borrow)

        # atilde_bc_borrow = astar_borrow - collateral * pbar * d_bc_borrow

        # === STEP 2(c): Interpolate onto a,d grid ===
        # interpolate onto a,d grid:
        c_borrow, d1_borrow, a_borrow, atilde_borrow = interp_to_base_grid( d_coh = d_bc_borrow,
                                                                            c_coh = c_bc_borrow)

        assert a_borrow.min()>-10**-8, 'error ' + np.array2string(a_borrow.min())

        # === STEP 3(a): Find solution for c',d' given tildea' = 0 (kink) ===
        # how to divide c',d' when a' is given

        # # NEED TO FIX KINKED SOLUTION

        d_bc_kink, c_bc_kink = solve_foc_x(astar = np.zeros([nE,nX]),
                                           interestrate_p = r_p ,
                                           dur_price = effective_dur_price_kink,
                                           dur_value_p = p_p * (1 - deltad))


        # === STEP 3(b): Interpolate onto a,d grid ===
        # interpolate onto a,d grid:
        c_kink, d1_kink, a_kink, atilde_kink = interp_to_base_grid( d_coh = d_bc_kink,
                                                                c_coh = c_bc_kink)


        # === STEP 4: Combine answers ===
        # realized borrowing / saving
        # print((atilde>0).sum(), (atilde_borrow<0).sum(), np.prod(atilde.shape) - (atilde>0).sum() - (atilde_borrow<0).sum())

        saving_region = (atilde>=0)
        borrowing_region = (atilde_borrow<=0)
        kink_region = (atilde<0) & (atilde_borrow>0)

        a_combine = np.empty([nE,nD,nA])
        a_combine[saving_region]    = a[saving_region]
        a_combine[kink_region]      = a_kink[kink_region]
        a_combine[borrowing_region] = a_borrow[borrowing_region]

        d1_combine = np.empty([nE,nD,nA])
        d1_combine[saving_region]    = d1[saving_region]
        d1_combine[kink_region]      = d1_kink[kink_region]
        d1_combine[borrowing_region] = d1_borrow[borrowing_region]

        c_combine = np.empty([nE,nD,nA])
        c_combine[saving_region]    = c[saving_region]
        c_combine[kink_region]      = c_kink[kink_region]
        c_combine[borrowing_region] = c_borrow[borrowing_region]

        a = a_combine
        d1 = d1_combine
        c = c_combine
    
    if spread>0:
        aborrow = np.zeros(a.shape)
        aborrow[borrowing_region] = a[borrowing_region]
    else:
        aborrow = np.zeros(a.shape)
     
    # === STEP 4: Update value function ===
    Vx = psi / c
    Va = (1 + interestrate0) * Vx
    Vd = current_dur_value * Vx

    V_points = np.concatenate((y_mesh.reshape([nN,1]), d1.reshape([nN,1]), a.reshape([nN,1])), axis=1)    
    V_p_interp = interpn((y, d1_grid, a_grid), V_p, V_points, method='linear', bounds_error=False, fill_value=None).reshape([nE,nD,nA])
    V = psi * np.log(c) + (1-psi) * np.log(d1) + beta * V_p_interp


    return Vx, Va, Vd, V, a, c, d1, aborrow



@het(exogenous='Pi', policy=['d1','a'], backward=['V','Vd','Va'], backward_init=household_init)
def household_localnofc(Va_p, Vd_p, V_p,  a_grid, d1_grid, x_grid, m_grid, s_grid, y, r, p, p_p, user_cost, op_cost, maint_cost, collateral, pbar, beta, eis, psi, deltad, xi, a_mesh, d1_mesh, y_mesh, y_mesh_fine, x_mesh_fine, s_mesh_fine, nE, nSfine, nA, nX, nD, dmax, spread):

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
    adj_dur_value = (p * (1 - deltad) - (1 + interestrate) * collateral * pbar)   

    util_params = {'eis': eis, 'xi': xi, 'psi': psi}
    
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
    # i_d, pi_d = interpolate.interpolate_coord(d1_grid, d_noadj.swapaxes(1, 2))

    Vd_p_inv_interp = interpolate.apply_coord(i_d, pi_d, (Vd_p ** -1).swapaxes(1, 2)).swapaxes(1, 2)
    Vd_p_inv_interp = interpolate.apply_coord(i_a, pi_a, Vd_p_inv_interp)
    Vd_p_interp = Vd_p_inv_interp ** -1

    V_p_interp = interpolate.apply_coord(i_d, pi_d, V_p.swapaxes(1, 2)).swapaxes(1, 2)
    V_p_interp = interpolate.apply_coord(i_a, pi_a, V_p_interp)


    Vd_noadj =  (1 - deltad_maint) * (f_ud(c_noadj, d_noadj, **util_params) - noadj_dur_value * Va_noadj + beta * Vd_p_interp)

    V_noadj = f_u(c_noadj, d_noadj, **util_params) + beta * V_p_interp
    
    # Va_noadj, Vd_noadj, V_noadj = update_value_noadj()
    
    # === ADJUSTER PROBLEM ===
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

    # === STEP 9: Update solutions and value function ===
    a = a_adj
    c = c_adj
    d1 = d_adj
    Va = Va_adj
    Vd = Vd_adj
    V = V_adj

    aborrow = np.zeros(a.shape)
    if spread>0:
        aborrow[borrowed] = - atilde_mesh[borrowed]

    
    assert np.isnan(V).sum()==0

    # check budget constraint
    coh_adj = (1 + interestrate) * a_mesh + adj_dur_value * d1_mesh + y_mesh
    assert np.max(np.abs(a_adj + effective_dur_price * d_adj + c_adj - coh_adj)) < 10**-8, 'budget_adj'

    
    return Va, Vd, V, a, c, d1, aborrow 