#%%
# import numpy as np
# import matplotlib.pyplot as plt

# these are the sequence space packges. they are in a separate directory for git management
import utils
from simple_block import simple
import jacobian as jac
import nonlinear


# ------------------------------------------------------------------------
# MODEL BLOCKS
# ------------------------------------------------------------------------

exogenous = dict()
unknowns = dict()
targets = dict()


@simple
def investment(k, iv, u, q, profitsk, lambdaa, rk, deltau, deltaudiff, adjC, adjCdiff, beta, delta, delta1, delta2, u_ss, kappa, gamma):
    # 0) Adjustment cost equation
    adjustcost = -adjC + kappa / 2 * (iv / iv(-1) - 1)**2
    adjustcostdiff = - adjCdiff + kappa * (iv / iv(-1) - 1)

    # 0) Depreciation
    depreciation = -deltau + delta + delta1 * (u - u_ss) + (delta2 / 2) * (u - u_ss)**2
    depreciationdiff = -deltaudiff + delta1 + delta2 * (u - u_ss)

    # 3) CONSUMPTION EULER EQUATION WITH TOBIN'S Q 
    eulerq = - lambdaa * q + beta * lambdaa(+1) * (rk(+1) * u(+1) + q(+1) * (1-deltau(+1)))

    # 4) INVESTMENT ADJUSTMENT COST
    investmentadj = (-lambdaa + lambdaa * q * (1 - adjC - (iv / iv(-1)) * adjCdiff) 
                + beta * lambdaa(+1) * q(+1) * (iv(+1) / iv)**2 * adjCdiff(+1) )

    # 5) OPTIMAL CAPITAL UTILIZATION
    optutil = -rk + q  * deltaudiff

    # 1) Capital equation of motion
    capitalacc = - k + (1-deltau)*k(-1) + iv*(1-adjC)

    # capital firm profits
    profk = - profitsk + rk * u * k(-1) - iv

    # input: rk, lambda
    # output: q, u, ivo, ko, iv, k

    return adjustcost, adjustcostdiff, depreciation, depreciationdiff, eulerq, investmentadj, optutil, capitalacc, profk

investment.exogenous = {'rk','lambdaa'}
investment.unknowns = {'q', 'u', 'iv', 'k', 'deltau', 'deltaudiff', 'adjC', 'adjCdiff', 'profitsk'}

# targets['investment'] = ['adjustcost', 'adjustcostdiff', 'depreciation', 'depreciationdiff', 
                        # 'eulerq', 'investmentadj', 'optutil', 'capitalacc', 'aggK', 'aggI']

    
@simple
def household_o(lambdaa, co, r, pi, ao, xo, do, dostar, go1, go2, w, ho, profitsk, to, profitso, pdur, tau, beta, sigma, gamma, deltad, psi, sigmad, varpi, xi, eta, thetad, incshock):
    # 2) SHADOW VALUE OF WEALTH
    if sigmad==sigma:
        if xi==1:
            complementarity = psi ** -psi * (1 - psi) ** (-(1 - psi)) * ( co ** (1 - psi) * do ** psi) ** (- 1 / sigma)
        else:
            complementarity = ((1 - psi)**(1 / xi) * co ** ((xi - 1) / xi) + psi ** (1 / xi) * do ** ((xi - 1) / xi)) ** (xi / (xi - 1) * (1 - 1 / sigma) - 1)

        shadowvalue = -lambdaa * tau + (1 - psi) ** (1/xi) * co **(-1/xi) * complementarity
        mud = psi ** (1/xi) * do **(-1/xi) * complementarity / tau
    else:
        shadowvalue = -lambdaa * tau + co **(-1/sigma)

    # 29) GROSS NOMINAL INTEREST RATE
    euler = -lambdaa + r * beta * lambdaa(+1) / pi(+1); 

    # RECURSIVE DURABLE CHOICE
    if sigmad==sigma:
        lomG1 = -go1 + mud * do ** (1 / varpi) + beta * thetad * (1 - deltad)  ** (1 - 1 / sigmad) * go1(+1)
    else:
        lomG1 = -go1 + psi + beta * thetad * (1 - deltad)  ** (1 - 1 / varpi) * go1(+1)
        
    lomG2 = -go2 + (pdur + eta) * tau * lambdaa - beta * (1 - deltad) * pdur(+1) * tau(+1) * lambdaa(+1) + beta * thetad * (1 - deltad) * go2(+1)


    # Durable RESET
    optd = -dostar + (go1 / go2) ** sigmad

    # DURABLE UPDATING
    updateD = -do + thetad * (1 - deltad) * do(-1) + (1 - thetad) * dostar

    # BUDGET CONSTRAINT
    budgeto = -ao/r + ao(-1)/pi + w * ho - tau * (co + eta * do + xo) - to + incshock + profitso + profitsk / (1 - gamma)

    # DURABLE EVOLUTION
    lomD = -do + (1 - deltad) * do(-1) + xo / pdur

    return shadowvalue, euler, budgeto, lomD, lomG1, lomG2, updateD, optd
    
household_o.exogenous = {'r','pi', 'b', 'to', 'ho', 'w', 'pdur', 'profitso', 'profitsk', 'tau', 'incshock'}
# household_o.unknowns = {'lambdaa', 'co', 'xo', 'do', 'ao'}
household_o.unknowns = {'lambdaa', 'co', 'xo', 'do', 'ao', 'go1', 'go2', 'dostar'}


@simple
def household_r(lambdaar, ar, cr, xr, dr, hr, r, pi, tr, w, tau, sigma, yr_ss, cr_ss, xr_ss, dr_ss, mpx, lag1, lag2, lag3, lag4, lag5, psi, pdur, deltad, eta, incshock):
    

    # 6) RULE-OF-THUMB CONSUMPTION
    RoTC = (- cr + ( cr_ss + 
           + (1-mpx) * (1 - lag1 - lag2 - lag3 - lag4 - lag5) * ( w*hr - tr - yr_ss + incshock) 
           + (1-mpx) * lag1 * ( w(-1)*hr(-1) - tr(-1) - yr_ss + incshock(-1)) * r(-1) / pi
           + (1-mpx) * lag2 * ( w(-2)*hr(-2) - tr(-2) - yr_ss + incshock(-2)) * r(-1) / pi * r(-2) / pi(-1)
           + (1-mpx) * lag3 * ( w(-3)*hr(-3) - tr(-3) - yr_ss + incshock(-3)) * r(-1) / pi * r(-2) / pi(-1) * r(-3) / pi(-2)
           + (1-mpx) * lag4 * ( w(-4)*hr(-4) - tr(-4) - yr_ss + incshock(-4)) * r(-1) / pi * r(-2) / pi(-1) * r(-3) / pi(-2) * r(-4) / pi(-3)
           + (1-mpx) * lag5 * ( w(-5)*hr(-5) - tr(-5) - yr_ss + incshock(-5)) * r(-1) / pi * r(-2) / pi(-1) * r(-3) / pi(-2) * r(-4) / pi(-3) * r(-5) / pi(-4)
            - eta * (dr - dr_ss) ) / tau )

    # 6) RULE-OF-THUMB Durable Expenditure
    RoTX = (- xr + ( xr_ss + 
           + mpx * (1 - lag1 - lag2 - lag3 - lag4 - lag5) * ( w*hr - tr - yr_ss + incshock) 
           + mpx * lag1 * ( w(-1)*hr(-1) - tr(-1) - yr_ss + incshock(-1)) * r(-1) / pi 
           + mpx * lag2 * ( w(-2)*hr(-2) - tr(-2) - yr_ss + incshock(-2)) * r(-1) / pi * r(-2) / pi(-1)
           + mpx * lag3 * ( w(-3)*hr(-3) - tr(-3) - yr_ss + incshock(-3)) * r(-1) / pi * r(-2) / pi(-1) * r(-3) / pi(-2)
           + mpx * lag4 * ( w(-4)*hr(-4) - tr(-4) - yr_ss + incshock(-4)) * r(-1) / pi * r(-2) / pi(-1) * r(-3) / pi(-2) * r(-4) / pi(-3)
           + mpx * lag5 * ( w(-5)*hr(-5) - tr(-5) - yr_ss + incshock(-5)) * r(-1) / pi * r(-2) / pi(-1) * r(-3) / pi(-2) * r(-4) / pi(-3) * r(-5) / pi(-4)
           ) / tau )

    # DURABLE EVOLUTION
    lomDr = -dr + (1 - deltad) * dr(-1) + xr / pdur           

    # Rule of thumb shadow value
    shadowvalue_r = - lambdaar * tau + cr**(-1/sigma)

    # BUDGET CONSTRAINT
    budgetr = -ar/r + ar(-1)/pi + w * hr - tau * (cr + eta * dr + xr) - tr

    return lomDr, RoTC, RoTX, shadowvalue_r, budgetr

household_r.exogenous = {'hr', 'tr','w', 'tau', 'r', 'pi', 'incshock'}
household_r.unknowns = {'dr', 'cr', 'xr', 'lambdaar', 'ar'}


@simple
def union(muw, lambdaa, lambdaar, h, ho, hr, hd, sw, f1, f2, w, wstar, pi, beta, gamma, nu, phi, thetaw, epsw, chiw):
    
    # lower part: (could swap hd or h)
    # input: cr, co, hr
    # output: f1, f2, wstar, hd, sw, w, h, ho, h, muw

    # 8) EQUALITY OF HOURS ACROSS TYPES
    equalH = -ho + hr

    # % 22) AGGREGATE HOURS
    aggH = -h + gamma * hr + (1 - gamma) * ho 

    # 7) MRS CONDITION
    MRS = - muw + w * (gamma / (nu * hr**phi)  * lambdaar + (1 - gamma) / (nu * ho**phi) * lambdaa)

    # 9)-11) SGU RECURSIVE EQS FOR OPTIMAL WAGE SETTING - TANK VERSION (substitutes gamma/cr + (1-gamma)/co for lambda in f1)

    focw1 = (-f1 + (epsw - 1) / epsw * wstar * (gamma  * lambdaar + (1-gamma) * lambdaa) * (w / wstar)**epsw * hd 
            + thetaw * beta * (pi(+1) / pi**chiw)**(epsw - 1) * (wstar(+1) / wstar)**(epsw - 1) * f1(+1) )
    focw2 = (-f2 + nu * h**phi * (w / wstar)**epsw * hd 
            + thetaw * beta * (pi(+1) / pi**chiw)**epsw * (wstar(+1) / wstar)**epsw * f2(+1) )
    focw = -f1 + f2

    # 12) - 13) DISTORTIONS DUE TO STICKY WAGES
    laborS = -h + sw * hd
    distortw = (-sw + (1 - thetaw) * (wstar / w)**(-epsw) 
                + thetaw * (w(-1 ) /w)**(-epsw) * (pi / pi(-1)**chiw)**epsw * sw(-1) )

    # 14) AGGREGATE WAGES
    wages = (-w**(1-epsw) + (1-thetaw) * wstar**(1-epsw) 
            + thetaw * w(-1)**(1-epsw) * (pi(-1)**chiw / pi)**(1 - epsw) )

    return MRS, equalH, aggH, focw1, focw2, focw, laborS, distortw, wages 
  
union.exogenous = {'lambdaar','lambdaa','hd'}
union.unknowns = {'f1', 'f2', 'wstar', 'hr', 'sw', 'w', 'h', 'ho', 'muw'}



@simple
def pricesetting(lambdaa, u, k, hd, w, rk, mc, x1, x2, p, pi, plvl, s, y, profits, profitso, alpha, beta, theta, eps, mup_ss, gamma):
    # 15) FOC FOR FIRM OPTIMAL CAPITAL-LABOR RATIO
    klratio = -u * k(-1) / hd + alpha / (1-alpha) * w / rk 

    # 16) REAL MARGINAL COST
    margC = -mc + rk**alpha * w**(1 - alpha) * alpha**(-alpha) * (1 - alpha)**(-1 + alpha) 

    # 17) - 19) URIBE RECURSIVE EQS FOR OPTIMAL PRICE SETTING (SGU 2005 SEEM DIFFERENT)
    focp1 = -x1 + mup_ss * p**(-eps) * y * mc + theta * beta * lambdaa(+1) / lambdaa * x1(+1)
    focp2 = -x2 + p**(1-eps) * y + theta * beta * lambdaa(+1) / lambdaa * p / p(+1) * x2(+1) / pi(+1)
    focp = -x1 + x2

    # 20) AGGREGATE INFLATION
    infl = -1 + (1 - theta) * p**(1 - eps) + theta * pi**(eps - 1)

    # price level
    pricelevel = - plvl + plvl(-1) * pi

    # 28) EQUATION OF MOTION FOR STICKY PRICE DISTORTION
    distortp = -s + (1 - theta) * p**(-eps) + theta * pi**eps * s(-1)

    # 27) PRODUCTION FUNCTION (S>=1 CAPTURES DISTORTION DUE TO STICKY PRICES)
    production = -s * y + (u * k(-1))**alpha * hd**(1-alpha) 

    # profit equations
    totprofits = - profits + y * (1 - mc)
    totprofitso = -profitso + profits / (1-gamma)

    # input: k, u, hd, w, lambda
    # output: rk, mc, x1, x2, p, pi, s, y

    return klratio, margC, focp1,focp2, focp, infl, distortp, production, pricelevel, totprofits, totprofitso

pricesetting.exogenous = {'k', 'u', 'y', 'w', 'lambdaa'}
pricesetting.unknowns = {'rk', 'mc', 'x1', 'x2', 'p', 'pi', 'plvl', 's', 'hd', 'profits', 'profitso'}



@simple
def monetary(r, pi, y, rr, rhor, r_ss, pi_ss, y_ss, phipi, phigap):

    # 30) INTEREST RATE RULE 
    intrule = -r + rhor * r(-1) + (1 - rhor) * (r_ss + phipi * (pi - pi_ss) + phigap * (y / y_ss - 1)) 

    # real rate
    realrate = -rr + r / pi(+1)

    # input: y, pi
    # output: r
    
    return intrule, realrate 

monetary.exogenous = {'pi', 'y'}
monetary.unknowns = {'r', 'rr'}



@simple
def fiscal(g, b, r, t, to, tr, pi, tau, c, x, tshock, g_ss, t_ss, to_ss, tr_ss, b_ss, exp_ss, phib, gamma):

    # 31) GOVERNMENT BUDGET CONSTRAINT
    govtbudget = -g + b/r - b(-1)/pi + t + (tau - 1) * (c + x)

    # 32) GOVERNMENT PURCHASES PROCESS
    spending = -g + g_ss

    # 33) - 34) TAX RULES - ASSUMES THEY DON'T BEGIN RAISING LUMP-SUM TAXES TO PAY OFF THE DEFICIT FOR MANY MONTHS
    taxrule = -t + t_ss + phib * (b(-12) - b_ss) - tshock * exp_ss
    taxdistr = -to + to_ss + tr - tr_ss # GLV implicit assumption (according to Uribe notes)

    # % 25) AGGREGATE TAXES
    aggT = -t + gamma * tr + (1-gamma) * to

    # input: r, pi
    # output: g, b, t, tr, to

    return govtbudget, spending, taxrule, taxdistr, aggT

fiscal.exogenous = {'r', 'tshock', 'tau', 'c', 'x'}
fiscal.unknowns = {'g', 'b', 't', 'tr', 'to'}



@simple
def marketclearing(c, co, cr, x, xo, xr, d, do, dr, nomc, nomx, nomexp, exp, plvl, pdur, realc, realx, iv, y, g, a, ao, ar, gamma, p_ss, eta):
    # % 21) AGGREGATE CONSUMPTION
    aggC = -c + gamma * cr + (1 - gamma) * co

    # % 21) AGGREGATE DURABLE EXPENDITURE
    aggX = -x + gamma * xr + (1 - gamma) * xo

    # % 21) AGGREGATE DURABLE STOCK
    aggD = -d + gamma * dr + (1 - gamma) * do

    REALX = - realx + x / pdur
    REALC = - realc + c + eta * d

    # % 26) RESOURCE CONSTRAINT
    aggY = -y + c + eta * d + iv + g + x

    #  TOTAL ASSETS
    aggA = -a + (1-gamma) * ao + gamma * ar

    # Nominal consumption
    nomC = -nomc + plvl * (c + eta * d)
    nomX = -nomx + plvl * x
    nomEXP = -nomexp + nomc + nomx
    EXP = -exp + c + eta * d + x / pdur
    # input: y, iv, g
    # output: c, cr

    return aggC, aggX, aggD, aggY, aggA, nomC, nomX, nomEXP, EXP, REALX, REALC

marketclearing.exogenous = {'dr', 'do', 'cr', 'co', 'xr', 'xo', 'iv', 'g', 'plvl', 'pdur'}
marketclearing.unknowns = {'d', 'c', 'x', 'y', 'a', 'nomc', 'nomx', 'nomexp', 'exp', 'realx', 'realc'}


@simple
def durablesupply(realx, pdur, zeta, realx_ss):
    # % Durable Supply elasticity
    durS = - pdur * realx_ss**(zeta) + realx**(zeta)

    # input: y, iv, g
    # output: c, cr

    return durS

durablesupply.exogenous = {'realx'}
durablesupply.unknowns = {'pdur'}

@simple
def twogoodsnondur(c1o, c2o, c1r, c2r, c1, c2, co, cr, cr_ss, share1, gamma):
    # optimizers consumption share
    ec1o = - c1o + share1 * co
    ec2o = - c2o + (1 - share1) * co

    ec1r = - c1r + (cr - cr_ss) + share1 * cr_ss
    ec2r = - c2r + (1 - share1) * cr_ss

    ec1 = - c1 + gamma * c1r + (1 - gamma) * c1o
    ec2 = - c2 + gamma * c2r + (1 - gamma) * c2o

    # input: y, iv, g
    # output: c, cr

    return ec1o, ec2o, ec1r, ec2r, ec1, ec2

twogoodsnondur.exogenous = {'ec1o', 'ec2o', 'ec1r', 'ec2r', 'ec1', 'ec2'}
twogoodsnondur.unknowns = {'c1o', 'c2o', 'c1r', 'c2r', 'c1', 'c2'}


for block in investment, household_r, household_o, union, pricesetting, monetary, fiscal, marketclearing, durablesupply:

        block.exogenous_list = list(block.exogenous)
        block.unknowns_list = list(block.unknowns)


# ------------------------------------------------------------------------
# TESTING BLOCKS
# ------------------------------------------------------------------------
#%%
# testing individual model blocks:
if __name__=='__main__':

    from nk_rebates_mpc_monthly_ss import nkrebate

    ss = nkrebate(thetad=0.1)
    T=50

    
    # checking individual blocks are correctly specified
    for block in investment, household_r, household_o, union, pricesetting, monetary, fiscal, marketclearing, durablesupply:
        

        assert len(block.unknowns_list) == len(block.output_list), 'elements in unknwons and targets do not match:' + str(block)

        

        G = jac.get_G(block_list=[block],
                    exogenous=block.exogenous_list,
                    unknowns=block.unknowns_list,
                    targets=block.output_list,
                    T=T, ss=ss)

        # print(G)

    # checking joing blocks are correctly specified
    G = dict()

    block_list = dict()
    block_list['PE'] = [household_o, household_r, marketclearing, fiscal]
    block_list['GE'] = [investment, household_o, union, pricesetting, monetary, fiscal, marketclearing, household_r, durablesupply]


    for key, blocks in block_list.items():
        allexogenous = set()
        allunknowns = set()
        targets = []
        for block in blocks:
            allexogenous = allexogenous | block.exogenous
            allunknowns = allunknowns | block.unknowns
            targets = targets + block.output_list


        exog = allexogenous - allunknowns
        unknowns = allunknowns - exog 
        unknowns = list( unknowns )
        exog = list( exog )

        G[key] = jac.get_G(block_list=blocks,
                        exogenous=exog,
                        unknowns=unknowns,
                        targets=targets,
                        T=T, ss=ss)

    # checking GE results clear asset market:
    assert (abs(G['GE']['a']['tshock']-G['GE']['b']['tshock'])).max() < 10**-6
    print('all checks passed')

    stop
    # comparing results with baseline Ramey model
    import pandas as pd
    variablenames = ['g', 'y', 'c', 'co', 'cr', 't', 'to', 'tr', 'h', 'ho', 'hr', 'w', 'iv', 'u', 'rk', 'r', 'pi', 'k', 'b'] #, , 'k', 'b']
    df = pd.read_csv('../output/testirfs.csv', names=variablenames)    

    print('Residuals: units in %')
    shock = 0.01
    for variable in variablenames:
        
        # print(G[variable]['tshock'][0:10,0],df[variable].iloc[0:10])
        if ss[variable]==0:
            steadydiv = ss['y']
        else:
            steadydiv = ss[variable]

        resid = ((G['GE'][variable]['tshock'][0:T,0] - df[variable][0:T]) / (steadydiv / shock + G['GE'][variable]['tshock'][0:T,0])).abs().max()  * 100
        print(variable + ': ' + str(resid))
        

