#%%
import numpy as np
import matplotlib.pyplot as plt


from sequence_jacobian import simple, solved, create_model


@solved(unknowns={'apc_inv': (-1e-6, 100000)}, targets=['apc_val'], solver="brentq")
def apc(sdf, apc_inv, beta, sigma, omega): 
    
    apc_val = -apc_inv + 1 + (beta ** sigma) * omega * sdf ** (1 - sigma) * apc_inv(+1)
    
    return apc_val

@solved(unknowns={'apd_inv': (-1e-6, 100000)}, targets=['apd_val'], solver="brentq")
def apd(sdf, Rd, apd_inv, beta, sigmad, omega, psi): 

    apd_val = -apd_inv + psi ** (sigmad) * Rd ** (1 - sigmad) + (beta ** sigmad) * omega * sdf ** (1 - sigmad) * apd_inv(+1)
    
    return apd_val    

@simple
def apcall(apc_inv, apd_inv, Rd, sigma, sigmad, psi, eta, C): 

    nondur_apc = (apc_inv + C ** (sigmad / sigma - 1) * apd_inv) ** -1
    dur_apc = psi ** (sigmad) * Rd ** ( - sigmad) * C ** (sigmad / sigma - 1) * nondur_apc
    oc_apc = eta * dur_apc
    nonduroc_apc = nondur_apc + oc_apc
    tot_apc = nondur_apc + dur_apc + oc_apc
    
    return nondur_apc, dur_apc, tot_apc, oc_apc, nonduroc_apc

@simple
def usercost(R, p, eta, deltad): 
    sdf = R(+1) ** -1
    Rd = p + eta - sdf * (1 - deltad) * p(+1)
    
    return Rd, sdf

@solved(unknowns={'Q': (-1000, 1000)}, targets=['asset_val'], solver="brentq")
def assets(Rd, R, D, X, B, Q, sdf, eta):
    Div = (Rd - eta) * D - X  + R * B(-1) - B     
    asset_val = -Q + sdf * ( Div(+1) + Q(+1) ) 
    Rs = (Q + Div) / Q(-1)
    
    return Div, Rs, asset_val

@solved(unknowns={'Qd': (-1000, 1000)}, targets=['acc'], solver="brentq")
def accumulation(Rs, Rd, Qd, Y, C, Q, D, T):
    acc = - Qd + Rs * Qd(-1) + Y - C - Rd * D - T
    asset_mkt = Qd - Q

    return acc, asset_mkt  

# @simple
# def accumulation(Rs, Rd, T, Y, C, Q, D):
#     Qd =  (Y - C - Rd * D - T) / (1 - Rs)
#     asset_mkt = Qd - Q

#     return Qd, asset_mkt        

@simple
def wealth_ss(sdf, Y, T, omega, Q, R):
    HW = (Y - T - (1 - omega) / omega * R * Q) /(1 - omega * sdf)
    FW = 1 / omega * R * Q
    
    return HW, FW      

@solved(unknowns={'HW': (0.01, 1000)}, targets=['human_wealth'], solver="brentq")
def wealth(sdf, HW, Y, Qd, Rs, T, omega, Q_ss, R_ss):
    human_wealth = - HW + Y - T - (1 - omega) / omega * R_ss * Q_ss + omega * sdf * HW(+1)
    FW = Rs * Qd(-1) + (1 - omega) / omega * R_ss * Q_ss

    return human_wealth, FW   

@simple
def budget_constr(C, D, HW, FW, Rd, nondur_apc, dur_apc): 
    
    budget_constr = - C + nondur_apc * (FW + HW )
    budget_constr_d = - D * Rd + dur_apc * (FW + HW )
    
    return budget_constr, budget_constr_d   

@simple
def euler(C, Q, Rs, nondur_apc, beta, sigma, omega, Q_ss): 
    
    euler_equ = -C + (beta * Rs(+1)) ** -sigma * ( C(+1) + (1 - omega) * nondur_apc(+1) * Rs(+1) / omega * (Q - Q_ss) )
    
    return euler_equ  

@simple
def euler_ss(C, R, beta, sigma): 
    
    euler_equ = -C + (beta * R(+1)) ** -sigma * C(+1)
    
    return euler_equ      

@simple
def optd(C, Rd, sigma, sigmad, psi): 
    D = psi ** (sigmad) * Rd ** (-sigmad) * C ** (sigmad / sigma)
    
    return D     

# @simple
# def optd_ss(C, Rd, sigma, sigmad, psi): 
#     D = psi ** (sigmad) * Rd ** (-sigmad) * C ** (sigmad / sigma)
    
#     return D   

@simple
def durable_price(p, X, X_ss, zeta): 
    dur_equ = -p * X_ss ** zeta + X ** zeta
    
    return dur_equ                
     
@simple
def aggregation(C, D, deltad, eta, p):
    X = (D - (1 - deltad) * D(-1)) * p
    expen = C + X + eta * D
    SX = X / (X + eta * D + C)
    
    return expen, X , SX  

@simple
def market_clearing(expen):
    Y = expen
    
    return Y   

@solved(unknowns={'B': (-10, 10)}, targets=['govtbudget'], solver="brentq")
def govt_budget(B, R, phib, T_ss, B_ss, T_shock):

    T = T_ss + phib * (B(-48) - B_ss) - T_shock

    govtbudget = B - R * B(-1) + T

    return govtbudget, T
        

blocks_ss = [usercost, apc, apd, apcall, assets, euler_ss, optd, aggregation, wealth_ss, accumulation, budget_constr, govt_budget, market_clearing]
blocks = [usercost, apc, apd, apcall, euler, optd, aggregation, wealth, accumulation, budget_constr, govt_budget, durable_price]
blocks_GE = blocks.copy()
blocks_GE.append(market_clearing)
blocks_GE.append(assets)

ppy = create_model(blocks_ss, name='PPY SS')

calibration = {'Y': 1, 'C': 1, 'R': 1.04**(1/12), 'R_ss': 1.04**(1/12),
                'omega': 0.975, 
                'psi': 0.00001, 'eta': 0.2/12, 'deltad': 0.015,
                 'p': 1, 'zeta': 0.2,
                'sigma': 0.2, 'sigmad': 1, 
                'B': 0, 'B_ss': 0, 'T': 0, 'T_ss': 0, 'T_shock': 0, 'phib':  0.01}

calibration['beta'] = calibration['R'] ** -1

unknowns_ss = {'C': 0.9, 'psi': 0.1}
targets_ss = {'Y': 1, 'SX': 0.04}
ss = ppy.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='broyden_custom')
# unknowns_ss = {'C': (0.8,1)}
# targets_ss = {'Y': 1}
# ss = ppy.solve_steady_state(calibration, unknowns_ss, targets_ss, solver='brentq')
# ss = ppy.steady_state(calibration)

# print(ss['tot_apc'])
# print(ss['nondur_apc'] + ss['oc_apc'])
# print(ss['dur_apc'])

assert abs((1-ss['deltad']) * ss['D'] / ss['R'] - ss['Q'])< 10 ** -8
assert abs(ss['budget_constr'])< 10 ** -8
assert abs(ss['euler_equ'])< 10 ** -8
assert abs(ss['asset_mkt'])< 10 ** -8
assert abs(ss['govtbudget'])< 10 ** -8


# assert abs(ss['goods_mkt'])< 10 ** -8

for var in ['C','Q','X', 'R']:
    ss[var + '_ss'] = ss[var]
ss['zeta'] = calibration['zeta']
ss['T_shock'] = calibration['T_shock']


ppy_micro = create_model(blocks, name='PPY micro')

exogenous = ['R', 'p', 'Y', 'T_shock']
unknowns = ['C']
targets = ['budget_constr']
T = 1000
G_micro = ppy_micro.solve_jacobian(ss, unknowns, targets, exogenous, T=T)

print(G_micro['expen']['Y'][0:3,0].sum())
print(G_micro['X']['Y'][0:3,0].sum())

print(G_micro['expen']['T_shock'][0:3,0].sum())
print(G_micro['X']['T_shock'][0:3,0].sum())
print(G_micro['X']['p'][0:6,6:T].sum() / (ss['X'] * 6))
print(G_micro['X']['R'][0:6,7].sum() / (ss['X'] * 6))

ppy = create_model(blocks_GE, name='PPY GE')
exogenous = ['R', 'T_shock']
unknowns = ['C','p']
targets = ['budget_constr', 'dur_equ']
T = 1000
G = ppy.solve_jacobian(ss, unknowns, targets, exogenous, T=T)

print(G['expen']['T_shock'][0:12,0].sum())
print(G['C']['T_shock'][0:12,0].sum())
print(G['X']['T_shock'][0:12,0].sum())
print(G['B']['T_shock'][-1,0])
plt.plot(G['B']['T_shock'][:,0])

# diffapc  = G['micro']['exp']['incshock'][0:3,0].sum() - apctarget
# diffelas = G['micro']['realx']['pdur'][0:lrhorizon,0:T['micro']].sum() / (ss['x'] * lrhorizon) - demandelasticty
# diffinter = G['X']['p'][0:6,6:T].sum() / (ss['X'] * 6) - 15

# def dag(cali, model = ppy):
    
#     cali['V_ss'] = (1 - cali['phi']) * cali['C'] / (cali['R'] - 1)

#     if cali['V_trans']==0:
#         unknowns_ss = {'beta': (1/cali['R'], 1/(cali['R'] * cali['omega'])),}
#         targets_ss = {'euler_equ': 0}

#         ss = model.solve_steady_state(cali, unknowns_ss, targets_ss, solver='brentq')

#     elif cali['V_trans']==1:
#         cali['beta'] = cali['R'] ** -1
#         ss = model.steady_state(cali)

#     exogenous = ['R']
#     unknowns = ['C']
#     targets = ['euler_equ']
#     T = 50

#     G = model.solve_jacobian(ss, unknowns, targets, exogenous, T=T)

#     return ss, G

# ss1, G1 = dag(calibration1)
# ss2, G2 = dag(calibration2)

# print('MPCs')
# print('With transfer: {}'.format(ss1['apc']))
# print('Without transfer: {}'.format(ss2['apc']))
# print('')
# print('Discount factor')
# print('With transfer: {}'.format(ss1['beta']))
# print('Without transfer: {}'.format(ss2['beta']))
# print('')
# print('Forward guidance')
# fg1 = G1['C']['R'][0,:10] * ss1['R']
# fg2 = G2['C']['R'][0,:10] * ss2['R']
# fig, ax = plt.subplots(1, 1, figsize=(6,4))
# ax.plot(np.linspace(0,10,10), fg1, label='Transfer', linewidth=2)
# ax.plot(np.linspace(0,10,10), fg2, label='No Transfer', linewidth=2)
# ax.legend(frameon=False)
# ax.axhline(0, color='gray', linestyle=':')
# ax.set_title('Forward Guidance')
# ax.set_xlabel('Horizon')
# ax.set_ylabel('Contemp Output')
# %%
