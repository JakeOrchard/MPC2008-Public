#%%
import numpy as np
# import matplotlib.pyplot as plt

# # these are the sequence space packges. they are in a separate directory for git management
# import sequence_jacobian_master.utils as utils
# import sequence_jacobian_master.jacobian as jac
# import sequence_jacobian_master.nonlinear as nonlinear
# from sequence_jacobian_master.simple_block import simple

# xc share: Table 2.4.5, PCE new vehicles / PCE 
# 233.2 / 9746.6 
# xc share: Table 2.4.5, PCE vehicle parts / PCE 
# 400.6 / 9746.6 = 0.0411
# alternative : PCE durable goods /  PCE
# 1188.0 / 9746.6  = 0.1219
# 2007 depreciation : all durables 879.8/4473.9
# vehicles 304.2 / 1330.7	
# maintenance (2017) vs 2016 stock: gas 319.9 + mv repair 143.4 + net insurance 62.0 / 2006 stock 1,316.5 / 12
# eta = (319.9+143.4+62)/1316.5 / 12= 0.03325
# main = (319.9+143.4+62) / 400.6  = 1.31

# ------------------------------------------------------------------------
# Steady state and parameter dictionary
# ------------------------------------------------------------------------

def nkrebate(delta   = 1-(1-.015)**(1/3),# monthly depreciation rate of private capital (annual = 6%, quarterly 1.5%)
            kappa    = 40,               # investment adj. cost, monthly value that approx IRF of 5.2 for quarterly from Leeper et al.
            phi      = 1,                # 1/phi is Frisch elasticity) GLV =  0.2
            alpha    = 0.36,             # capital share, GLV = 1/3
            mup      = 1.2,              # steady state price markup
            muw      = 1.2,              # wage markup
            theta    = 0.9167,           # price stickiness, avg. duration of a price = 1/(1-theta). monthly equiv to theta = 0.75 at quarterly, 0 implies no stickiness, GLV = 0.75
            thetaw   = 0.9167,           # wage stickiness, monthly equiv to theta = 0.75 at quarterly, 0 implies no stickiness
            chiw     = 0,                # wage indexing, , 0=no, 1=full
            sigma    = 0.5,              # IES
            gamma    = 0,                # rule-of-thumb share, GLV = 0.5 
            mpx      = 0.83,             # rule-of-thumb marginal expenditure on durables
            lag1     = 1/3,                # rule-of-thumb mpc on transfer last month
            lag2     = 1/3,                # rule-of-thumb mpc on transfer two months ago
            lag3     = 0,                # rule-of-thumb mpc on transfer three months ago
            lag4     = 0,                # rule-of-thumb mpc on transfer four months ago
            lag5     = 0,                # rule-of-thumb mpc on transfer five months ago
            omega    = 1,                # survival probability
            r        = 0.99**(-1/3),     # monthly gross interest rate
            xc       = 0.0411,           # ratio of durable stock to nondurable expenditure
            deltad   = 1.2**(1/12)-1,    # depreciation rate of durable consumption stock
            sigmad   = 1,                # exponent on durable service flow in utility function
            xi       = 1,                # elasticity of substitution between nondurable and durable goods
            varpi    = 1,                # elasticity of substitution between durable goods varieties
            thetad   = 0,                # durable stickiness
            maintx   = 1.31,            # ratio of maintenance / operating cost to car expenditures
            eta      = 0.22 / 12,
            zeta     = 0,                # inverse durable supply elasticity
            phib     = 0.027,              # fiscal rule parameters, our quarterly = 0.1, GLV quarterly = 0.33 
            phig     = 0,                # fiscal rule parameters, lower phig means less increase in tax, GLV = 0.1 
            phipi    = 1.5,              # Taylor rule parameter on inflation
            phigap   = 1/12,             # *** Taylor rule parameter on output gap, 1 for annual in balanced-approach rule, divide by 12 for monthly
            rhor     = 0.85**(1/3),      # Inertia in Taylor rule, 0 = no inertia, quarterly rhor = monthly rhor^3, .85 from Fed website
            pi       = 1,                # steady state inflation
            y        = 1,                # steady state level of output
            gyfrac   = 0.2,              # Steady-state G/Y, set to 0, not important in this project
            rhog     = 0,                # AR(1) on government spending
            ubar     = 1,                # steady-state utilization
            eutil    = 2,                # marginal depreciation elasticity to utilization
            share1   = 0.528,            # share of type 1 good in nondurable expenditure
            ):
    """Solve steady state of NK Rebate model.
    Returns
    -------
    ss : dict, steady state values
    """

    # Implied parameters
    eps      = mup / (mup - 1) # elasticity of subsitution of goods varieties
    epsw     = muw / (muw - 1) # elasticity of substitution of labor varieties
    delta1   = r - 1 + delta   # utiliz. cost linear term, so u = ubar in SS
    delta2   = eutil*delta1               # utiliz. cost quadratic term, SGU 2.02 based on ACEL, Leeper et al. approx. 2

    # solve for aggregates analytically
    gy = gyfrac        # Steady-state g/y ratio
    s = 1
    p = 1
    pdur = 1
    beta = 1 / r
    rr = r / pi
    rk = r -1+delta
    mc = 1/mup
    sw = 1
    u = ubar
    q = 1
    b = 0
    a = 0
    tau = 1
    
    xd = 1 - (1 - deltad)
    ky = alpha / ( rk * mup)
    iy = delta * ky
    gy = gyfrac
    # cy = ( 1 - iy - gy ) / (1 + xc * (1 + maintx))
    # eta = 0.04
    # wealth = (1 - omega) * (a + (1 - deltad)
    cy = ( 1 - iy - gy ) / (1 + xc * (1 + eta / xd))
    dc = xc * xd
    kh = ky ** (1 / (1-alpha))
    w = kh * rk * (1-alpha) / alpha

    k = ky * y
    iv = iy * y
    g = gy * y
    c = cy * y
    x = xc * c
    d = x / xd
    h = k / kh
    
    # implied disutility of labor
    xo = x
    xr = x
    do = d
    dr = d
    co = c
    cr = c

    # calvo parameters
    dostarratio = (1 - thetad * (1 - deltad)) / (1 - thetad)
    dostar = dostarratio * do
    drstar = dostar

    # marginal utility of consumption
    if do==0:
        eta = 0
        psi = 0
            
    user_cost = p + eta - p * (1 - deltad) / r
    rd = user_cost * (1 - beta * thetad * (1 - deltad)**(1 - 1 / varpi) ) / (1 - beta * thetad * (1 - deltad)) * dostarratio ** (1 / varpi)


    muc = co ** (-1 / sigma)

    # recursive formulation of durable optimization
    go2 = muc * user_cost / (1 - beta * thetad * (1 - deltad))

    go1 = dostar ** (1 / sigmad) * go2    
    psi = go1 * (1 - beta * thetad * (1 - deltad)  ** (1 - 1 / sigmad)) 

    if psi>0:
        mud = psi * do ** (-1/sigmad)
    else:
        mud = 0
        
    if sigma==sigmad and do>0:
        psi = (do / co) / (do / co + rd**(-xi))

        if xi!=1:
            complementarity = ((1 - psi)**(1 / xi) * co ** ((xi - 1) / xi) + psi ** (1 / xi) * do ** ((xi - 1) / xi)) ** (xi / (xi - 1) * (sigma - 1) / sigma - 1)
        else:
            complementarity = psi ** -psi * (1 - psi) ** (-(1 - psi)) * ( co ** (1 - psi) * do ** psi) ** (- 1 / sigma)

        muc = (1 - psi)**(1 / xi) * co ** (- 1 / xi) * complementarity
        mud =      psi **(1 / xi) * do ** (- 1 / xi) * complementarity
        
        lambdatimesrd = mud

        # recursive formulation of durable optimization
        go2 = muc * user_cost / (1 - beta * thetad * (1 - deltad))
        go1 = muc * rd * do ** (1 / varpi) / (1 - beta * thetad * (1 - deltad)**(1 - 1 / varpi))

        assert np.isclose(dostar, (go1/go2)**varpi ), "complementarity equations not correct"
            

    # marginal utility of durable good
    lambdaa = muc
    


    mucr = muc
    mudr = mud
    lambdaar = mucr

    lambdatilde = (1 - gamma) * lambdaa + gamma * lambdaar

    nu = lambdatilde * w / (muw * h**phi)

    
    hd = h
    wstar = w
    ho = h
    hr = h
    ko = k / (1-gamma)
    ivo = iv / (1 - gamma)
    ao = a / (1-gamma)
    ar = 0
    deltau = delta
    deltaudiff = delta1
    adjC = 0
    adjCdiff = 0
    tr = w * hr - cr - xr - eta * dr
    to = (g - gamma*tr) / (1-gamma)
    t = (1-gamma) * to + gamma * tr
    yr = w * hr - tr
    tshock = 0
    incshock = 0
    x1 = y / (1 - theta * beta)
    x2 = x1
    f1 = w * lambdaa * hd * muw**(-1) / (1 - thetaw * beta)
    f2 = f1
    plvl = 1
    nomc = plvl * (c + eta * d)
    nomx = plvl * x
    nomexp = nomc + nomx
    exp = c + eta * d + x / p
    realx = x / p
    realc = c + eta * d
    profitsk = rk * u * k - iv
    profits  = (1 - 1/mup) * y
    profitso = profits / (1-gamma)
    
    # if there are no durables set mpx=0:
    if xc==0:
        mpx = 0 
        zeta = 0
        sigmad = 1
        eta = 0
        # print('No Durables in Steady State. Setting MPX = 0 and supply elasticity = 0.')

    # variables for 2001 calibration
    c1r = share1 * cr
    c2r = (1 - share1) * cr
    c1o = share1 * co
    c2o = (1 - share1) * co
    c1 = gamma * c1r + (1 - gamma) * c1o
    c2 = gamma * c2r + (1 - gamma) * c2o

    res = dict()
    res['budgetr'] = - cr - eta * dr - xr + w * hr - tr
    res['budgeto'] = - co - eta * do - xo + w * ho - to + profitsk / (1-gamma) + profitso
    res['kh'] = -kh + k / h
    res['capital accum'] = -k + (1- delta) * k + iv * (1 - kappa/2 * (iv/iv -1)**2 )
    res['shadow value'] = lambdaa - muc
    res['shadow value r'] = lambdaar - mucr
    # res['relative durable'] = lambdaa * (1 + eta) * do ** (1 / sigmad) - beta * (1 - deltad) * (1 - fc) * lambdaa * do ** (1 / sigmad) - psi 
    # res['relative durable'] = go2 * do ** (1 / sigmad) - psi
    res['calvo durable'] = do - (1 - thetad) * dostar - thetad * (1 - deltad) * do
    if sigmad==sigma:
        res['optimal reset'] = dostar - (go1 / go2) ** varpi
        res['relative durable'] = lambdaa * p * dostar ** (1 / varpi) + lambdaa * eta / (1 - beta * thetad * (1 - deltad)) * dostar ** (1 / varpi) - beta * (1 - thetad) * (1 - deltad) / (1 - beta * thetad * (1 - deltad)) * lambdaa * dostar ** (1 / varpi) - lambdatimesrd * do ** (1 / varpi) / (1 - beta * thetad * (1 - deltad) ** (1 - 1 / varpi))
        res['calvo numerator'] = go1 - lambdatimesrd * do ** (1 / varpi) - beta * thetad * (1 - deltad)  ** (1 - 1 / varpi) * go1
    else:
        res['optimal reset'] = dostar - (go1 / go2) ** sigmad
        res['relative durable'] = lambdaa * p * dostar ** (1 / sigmad) + lambdaa * eta / (1 - beta * thetad * (1 - deltad)) * dostar ** (1 / sigmad) - beta * (1 - thetad) * (1 - deltad) / (1 - beta * thetad * (1 - deltad)) * lambdaa * dostar ** (1 / sigmad) - psi / (1 - beta * thetad * (1 - deltad) ** (1 - 1 / sigmad))
        res['calvo numerator'] = go1 - psi - beta * thetad * (1 - deltad)  ** (1 - 1 / sigmad) * go1
    res['calvo denominator'] = go2 - lambdaa * (1 + eta - beta * (1 - deltad)) - beta * thetad * (1 - deltad) * go2
    res['durable acc'] = do  - (1 - deltad) * do - xo
    # print(lambdaa)
    # print(mud)
    # print(muc)
    # res['relative durable r'] = lambdaar * (1 + eta) * dr ** (1 / varphi) - beta * (1 - deltad) * lambdaar * dr ** (1 / varphi) - psi ** (1 / varphi) * complementr
    # res['relative durable r'] = lambdaar * (1 + eta) * dr ** (1 / sigmad) - beta * (1 - deltad) * (1 - fc) * lambdaar * dr ** (1 / sigmad) - psi
    res['durable acc r'] = dr - (1 - deltad) * dr - xr
    res['euler equation'] = lambdaa * q - beta * lambdaa * (rk * u + q * (1 - delta))
    res['investment adj C'] = - lambdaa + ( lambdaa * q  )
    res['optimal capital acc'] = rk - q * (delta1 + delta2*(u-1))
    res['HTM consumers'] = - cr - eta * dr - xr + w * hr - tr
    res['MRS'] = - muw + w * (gamma / (nu * hr ** phi) * lambdaar + (1-gamma) / (nu * ho ** phi) * lambdaa)
    res['hours equality'] = ho - hr
    res['f1'] = ( -f1 + (epsw-1)/epsw * wstar * (gamma * lambdaar + (1-gamma) * lambdaa) * (w/wstar)**epsw * hd 
                + thetaw * beta * pi ** ((1-chiw)*(epsw-1)) * f1 )
    res['f2'] = (-f2 + nu * h**phi * (w/wstar) ** epsw * hd + thetaw * beta * pi**((1-chiw)*epsw) * f2)
    res['ff'] = f1 - f2
    res['dist1'] = -h + sw * hd
    res['dist2'] = -sw + (1-thetaw) * (wstar/w)**(-epsw) + thetaw* pi**((1-chiw)*epsw) * sw
    res['wages'] = -w**(1-epsw) + (1-thetaw) * wstar**(1-epsw) + thetaw * w**(1-epsw) * pi**((chiw-1)*(1-epsw))
    res['kl ratio'] = - u * k / hd + alpha/(1-alpha) * w / rk
    res['real MC'] = -mc + rk**alpha * w**(1-alpha) * alpha**(-alpha) * (1-alpha)**(-1+alpha)
    res['x1'] = -x1 + mup * p**(-eps) * y * mc + theta * beta * x1
    res['x1'] = -x2 + p**(1-eps) * y + theta * beta * x2 / pi
    res['xx'] = -x1 + x2
    res['inflation'] = -1 + (1-theta) * p**(1-eps) + theta * pi**(eps-1)
    res['consumption'] = -c + gamma*cr + (1-gamma)*co
    res['durale exp'] = -x + gamma*xr + (1-gamma)*xo
    res['durale stock'] = -d + gamma*dr + (1-gamma)*do
    res['hours'] = -h + gamma*hr + (1-gamma)*ho
    res['investment'] = -ivo + iv / (1-gamma)
    res['capital'] = -ko + k / (1-gamma)
    res['taxes'] = -t + gamma*tr + (1-gamma)*to
    res['resources'] = -y + c + eta * d + iv + g + x
    res['production'] = -s * y + (u * k) ** alpha * hd ** (1-alpha)
    res['sticky P distortion'] = -s + (1-theta) * p**(-eps) + theta * pi**eps * s
    res['nominal rate'] = -lambdaa + r * beta * lambdaa / pi
    res['sticky P distortion'] = -s + (1-theta) * p**(-eps) + theta * pi**eps * s
    res['interest rate rule'] = -r + rhor * r + (1 - rhor) * r 
    res['govt budget constraint'] = -g + b/r - b/pi + t
    res['tax rules'] = -t + t + phib*(b - b) - tshock

    


    


    steady = {                                                      # VARIABLES
            'co': co, 'cr': cr, 'c': c,                             # consumption
            'xo': xo, 'xr': xr, 'x': x,                             # durable expenditure 
            'nomc': nomc, 'nomx': nomx, 'nomexp': nomexp,           # nominal expenditure
            'exp': exp, 'realx': realx, 'realc': realc,             # real expenditure
            'do': do, 'dr': dr, 'd': d, 'pdur': pdur,               # durable stock
            'dostar': dostar, 'go1': go1, 'go2': go2,               # durable updating
            'ho': ho, 'hr': hr, 'h': h, 'hd': hd,                   # hours
            'k': k, 'iv': iv, 'q': q, 'u': u, 'profitsk': profitsk, # investment and capital stock, utilization
            'deltau': deltau, 'deltaudiff': deltaudiff,             # depreciation function 
            'adjC': adjC, 'adjCdiff': adjCdiff,                     # adjustment cost function
            'y': y, 'yr': yr,                                       # output
            'a': a, 'ao': ao, 'ar': ar,                             # assets
            'lambdaa': lambdaa, 'lambdaar': lambdaar, 'mud': mud,   # Lagrange multiplier on HH budget constraint, 
            'w': w, 'rk': rk,                                       # wages, rental rate of capital
            'p': p, 'mc': mc, 'muw': muw,                           # relative price, marginal cost, wage markup 
            's': s, 'sw': sw, 'wstar': wstar,                       # relative price, marginal cost, wage markup 
            'x1': x1, 'x2': x2, 'f1': f1, 'f2': f2,                 # auxiliary terms for price and wage stickiness
            'profits': profits, 'profitso': profitso,               # profits
            'r': r, 'pi': pi, 'rr': rr, 'plvl': plvl,               # gross nominal interest rate,inflation
            't': t, 'tr': tr, 'to': to, 'b': b, 'tshock': tshock,   # total taxes, taxes on ROT HH, taxes on optimizing HH, govt debt  
            'g': g, 'tau': tau, 'incshock': incshock,               # government purchases and sales taxes
            'c1r': c1r, 'c2r': c2r, 'c1o': c1o, 'c2o': c2o,         # consumption by type
            'c1': c1, 'c2': c2,                                     # consumption by type
                                                                    # PARAMETERS: STEADY STATE
            'co_ss': co, 'cr_ss': cr, 'c_ss': c,                    # consumption
            'xo_ss': xo, 'xr_ss': xr, 'x_ss': x,                    # durable expenditure 
            'nomc_ss': nomc, 'nomx_ss': nomx, 'nomexp_ss': nomexp,  # nominal expenditure
            'exp_ss': exp, 'realx_ss': realx,                       # real expenditure
            'do_ss': do, 'dr_ss': dr, 'd_ss': d, 'pdur_ss': pdur,   # durable stock
            'ho_ss': ho, 'hr_ss': hr, 'h_ss': h, 'hd_ss': hd,       # hours
            'k_ss': k, 'iv_ss': iv, 'q_ss': q, 'u_ss': u, 'profitsk_ss': profitsk, # investment and capital stock, utilization
            'deltau_ss': deltau, 'deltaudiff_ss': deltaudiff,       # depreciation function 
            'adjC_ss': adjC, 'adjCdiff_ss': adjCdiff,               # adjustment cost function
            'y_ss': y, 'yr_ss': yr,                                 # output
            'a_ss': a, 'ao_ss': ao,                                 # assets
            'lambdaa_ss': lambdaa, 'lambdaar_ss': lambdaar, 'mud_ss': mud, # Lagrange multiplier on HH budget constraint,
            'w_ss': w, 'rk_ss': rk,                                 # wages, rental rate of capital
            'p_ss': p, 'mc_ss': mc, 'muw_ss': muw,                  # relative price, marginal cost, wage markup 
            's_ss': s, 'sw_ss': sw, 'wstar_ss': wstar,              # relative price, marginal cost, wage markup 
            'x1_ss': x1, 'x2_ss': x2, 'f1_ss': f1, 'f2_ss': f2,     # auxiliary terms for price and wage stickiness
            'profits_ss': profits, 'profitso_ss': profitso,         # profits
            'r_ss': r, 'pi_ss': pi, 'rr_ss': rr, 'plvl_ss': plvl,   # gross nominal interest rate,inflation
            't_ss': t, 'tr_ss': tr, 'to_ss': to, 'b_ss': b , 'tshock_ss': tshock, # total taxes, taxes on ROT HH, taxes on optimizing HH, govt debt  
            'g_ss': g, 'tau': tau,                                  # government purchases and sales taxes
            'c1r_ss': c1r, 'c2r_ss': c2r, 'c1o_ss': c1o, 'c2o_ss': c2o, # consumption by type
            'c1_ss': c1, 'c2_ss': c2,                               # consumption by type
              
                                                                    # PARAMETERS: TECHNOLOGY, PREFERENCES

            'beta': beta,                                           # discount factor
            'sigma': sigma,                                         # IES
            'psi': psi,                                             # weight on durables in utility
            'sigmad': sigmad, 
            'varpi': varpi,
            'xi': xi,                                                                            # curvature of durables in utility
            'thetad': thetad,                                       # calvo probability for durable choice
            'deltad': deltad,                                       # durable depreication rate
            'eta': eta,                                             # durable operating cost
            'phi': phi,                                             # utility parameter (1/phi is Frisch elasticity)
            'nu': nu,                                               # weight on disutility of labor
            'gamma': gamma,                                         # rule-of-thumb share
            'mpx': mpx,                                             # rule-of-thumb marginal expenditure on durables
            'lag1': lag1,                                            # rule-of-thumb mpc on transfer last month
            'lag2': lag2,                                            # rule-of-thumb mpc on transfer two months ago
            'lag3': lag3,                                            # rule-of-thumb mpc on transfer three months ago
            'lag4': lag4,                                            # rule-of-thumb mpc on transfer four months ago
            'lag5': lag5,                                            # rule-of-thumb mpc on transfer five months ago
            'delta': delta,                                         # depreciation rate on private capital                                 
            'kappa': kappa,                                         # investment adjustment cost
            'delta1': delta1,                                       # parameter of utilization cost
            'delta2': delta2,                                       # parameter of utilization cost
            'ubar': ubar,                                           # steady-state capital utilization
            'alpha': alpha,                                         # capital share
            'zeta': zeta,                                           # durable supply elasticity
            'theta': theta,                                         # price stickiness
            'thetaw': thetaw,                                       # wage stickiness
            'chiw': chiw,                                           # degree of wage indexing
            'phipi': phipi,
            'phigap': phigap,                                       # Taylor rule parameters
            'phig': phig,
            'phib': phib,                                           # fiscal financing parameters
            'eps': eps,                                             # elasticity of substitution of goods varieties
            'epsw': epsw,                                           # elasticity of substitution of labor varieties
            # 'muw_ss': muw,                                          # steady-state wage markup
            'mup_ss': mup,                                          # steady-state price markup
            'gyfrac': gyfrac,                                       # G/Y fraction
            'rhog': rhog,                                           # AR(1) on government spending
            'rhor': rhor,  
            'share1': share1,                                           # Inertia in the Taylor Rule
           }

    for key, value in res.items():
        assert np.abs(value) < 1E-5, 'Error in SS equation for ' + key + ', value = ' + str(value)          

    return steady

# test that number of elements is as expected
if __name__=='__main__':
    # ss = nkrebate(zeta=0.2, gamma=0.5, varphi=0.8, sigma=0.5, sigmad=0.5, elasdv=0.5, maintx=1.3, fc=0.05)

    # ss = nkrebate(zeta=0.2, gamma=0.5, varphi=0.8, sigma=0.5, sigmad=0.5, elasdv=0.5, maintx=1.3, thetad=0, fc=0)

    # ss = nkrebate(zeta=0.2, gamma=0.5, varphi=0.8, sigma=0.5, sigmad=0.5, elasdv=0.5, maintx=1.3, thetad=0.5, fc=0)

    # ss = nkrebate(zeta=0.2, gamma=0.5, varphi=0.8, sigma=0.5, sigmad=0.5, elasdv=0.25, maintx=1.3, thetad=0.5, fc=0)

    ss = nkrebate(zeta=0.2, gamma=0.5, sigmad=0.8, sigma=0.5, xc=0, thetad=0)

    ss = nkrebate(zeta=0.2, gamma=0.5, sigmad=0.8, sigma=0.5, xc=0, thetad=0.5)
    # assert len(ss)== 56*2 + 30
    ss = nkrebate(zeta=0.2, gamma=0.5, sigmad=0.8, sigma=0.5, thetad=0.5)

    ss = nkrebate(zeta=0.2, gamma=0.5, sigmad=0.5, sigma=0.5, thetad=0.5)

    ss = nkrebate(zeta=0.2, gamma=0.5, sigmad=0.5, sigma=0.5, thetad=0.5, xi=0.5)

    # test2 = nkrebate(sigma=2)

    # test3 = nkrebate(gamma=0)

# %%
