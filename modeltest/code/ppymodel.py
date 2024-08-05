#%%
import numpy as np
import matplotlib.pyplot as plt


from sequence_jacobian import simple, solved, create_model


@solved(unknowns={'mpc_inv': (1e-10, 1000)}, targets=['mpc_val'], solver="brentq")
def mpccalc(R, mpc_inv, beta, sigma, omega): 
    sdf = R(+1) ** -1
    mpc_val = -mpc_inv + 1 + (beta ** sigma) * omega * (sdf ** (1 - sigma)) * mpc_inv(+1)
    mpc = mpc_inv ** -1
    
    return mpc_val, sdf, mpc


@solved(unknowns={'Q': (-0.00001, 20000)}, targets=['asset_val'], solver="brentq")
def assets(sdf, phi, Y, Q):
    asset_val = -Q + sdf * ( (1 - phi) * Y(+1) + Q(+1) ) 
    Div = (1 - phi) * Y
    Rs = (Q + Div) / Q(-1)
    
    return asset_val, Div, Rs

@simple
def euler(C, Qd, R, Rs, mpc, beta, sigma, omega, Q_ss, V_trans, R_ss): 
    euler_equ = -C + (beta * R(+1)) ** -sigma *  C(+1) + (beta * R(+1)) ** -sigma * (1 - omega) * mpc(+1)  / omega * (Rs(+1) * Qd - R_ss * V_trans * Q_ss) 
    
    return euler_equ    
     
@simple
def market_clearing(Y): 
    C = Y
    
    return C 

@solved(unknowns={'Qd': (-0.00001, 20000)}, targets=['acc'], solver="brentq")
def accumulation(Rs, Qd, Y, C, Q, phi):
    acc = - Qd + Rs * Qd(-1) + phi * Y - C
    asset_mkt = Qd - Q

    return acc, asset_mkt

@solved(unknowns={'HW': (-100, 100000)}, targets=['human_wealth'], solver="brentq")
def wealth(sdf, HW, Y, Rs, Qd, omega, phi, Q_ss, R_ss, V_trans):
    Transfer = V_trans * (1 - omega) / omega * R_ss * Q_ss
    
    human_wealth = - HW + phi * Y - Transfer + omega * sdf * HW(+1)
    FW = Rs * Qd(-1)

    return human_wealth, FW, Transfer

@simple
def budget_constr(C, HW, FW, mpc, Transfer, omega, R_ss, V_trans, Q_ss): 
    
    budget_constr = - C + mpc * (FW + HW  + Transfer )
    impl_mpc = C / (FW + HW  + Transfer )
    diff_c = mpc  * (omega ** -1) * (FW - R_ss * V_trans * Q_ss) 
    
    return budget_constr, impl_mpc, diff_c 

blocks_micro = [mpccalc, assets, euler, wealth, budget_constr, accumulation]
blocks_GE = [mpccalc, assets, euler, market_clearing, wealth, budget_constr, accumulation]

ppy = create_model(blocks_GE, name='PPY')
ppy_micro = create_model(blocks_micro, name='PPY micro')

calibration1 = {'C': 1., 'Y': 1., 'R': 1.01, 'R_ss': 1.01,
                'omega': 0.975, 'phi': 0.9,
                'sigma': 0.5, 'V_trans': 1}

calibration2 = calibration1.copy()
calibration2['V_trans'] = 0

def dag(cali, model = ppy, model_micro = ppy_micro):
    
    cali['Q_ss'] = (1 - cali['phi']) * cali['C'] / (cali['R'] - 1)

    if cali['V_trans']<1:
        unknowns_ss = {'beta': (1/cali['R']-0.01, (1/(cali['R'] * cali['omega']))**(1/cali['sigma'])+0.01),}
        targets_ss = {'euler_equ': 0}

        ss = model.solve_steady_state(cali, unknowns_ss, targets_ss, solver='brentq')

    elif cali['V_trans']==1:
        cali['beta'] = cali['R'] ** -1
        ss = model.steady_state(cali)

    exogenous = ['R']
    unknowns = ['Y']
    targets = ['budget_constr']
    T = 50

    G = model.solve_jacobian(ss, unknowns, targets, exogenous, T=T)

    exogenous = ['Y', 'R']
    unknowns = ['C']
    targets = ['budget_constr']
    T = 50
    G_micro = model_micro.solve_jacobian(ss, unknowns, targets, exogenous, T=T)
    

    return ss, G, G_micro

ss1, G1, G1_micro = dag(calibration1)
ss2, G2, G2_micro = dag(calibration2)

assert abs(ss1['budget_constr'])<10**-8
assert abs(ss2['budget_constr'])<10**-8
assert abs(ss1['acc'])<10**-8
assert abs(ss2['acc'])<10**-8

assert abs(G1['budget_constr']['R']).max()<10**-8
assert abs(G2['budget_constr']['R']).max()<10**-8

assert abs(G1['euler_equ']['R']).max()<10**-8
assert abs(G2['euler_equ']['R']).max()<10**-8

assert abs(G1['acc']['R']).max()<10**-8
assert abs(G2['acc']['R']).max()<10**-8



print('MPCs')
print('With transfer: {}'.format(ss1['mpc']))
print('Without transfer: {}'.format(ss2['mpc']))
print('')
print('Discount factor')
print('With transfer: {}'.format(ss1['beta']))
print('Without transfer: {}'.format(ss2['beta']))
print('')
print('Forward guidance')
fg1 = G1['C']['R'][0,1:11] * ss1['R']
fg2 = G2['C']['R'][0,1:11] * ss2['R']
fig, ax = plt.subplots(1, 1, figsize=(6,4))
ax.plot(np.linspace(0,10,10), fg1, label='Transfer', linewidth=2)
ax.plot(np.linspace(0,10,10), fg2, label='No Transfer', linewidth=2)
ax.legend(frameon=False)
ax.axhline(0, color='gray', linestyle=':')
ax.set_title('Forward Guidance')
ax.set_xlabel('Horizon')
ax.set_ylabel('Contemp Output')

# %%
