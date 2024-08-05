#%%
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from cycler import cycler

plt.rc('font', size=16) #controls default text size

# ------------------------------------------------------------------------
# Tables and Figure Settings
# ------------------------------------------------------------------------

parasetnames = ['baseline','durableprice','lesselastic','nondurablesonly','sameexponent','sameexponentdurableprice','sameexponentdurablepricehigh','sameexponentdurablepricelow','inelasticdurables'] # 'gefeedback',,'inelasticsupply','strongmp' #['inelasticdurables']#
parasettitles = ['Baseline Model','Model with Less Elastic Durable Supply','Model with Less Elastic Durable Demand','Model without Durable Goods','Same Exponent Model','Same Exponent Less Elastic Durable Supply','Same Exponent Less Elastic Durable Supply High Substitutability','Same Exponent Less Elastic Durable Supply Low Substitutability','Model without Inelastic Durable Supply']

# parasetnames = ['roffset']
# parasettitles = ['Model with Constant Output Gap']



plot_names = ['PCE', 'Motor vehicles', 'Nondurable goods', 'Relative Durable Price', 'Federal Funds Rate'] #,

table_dict = {'PCE': 'exp', 'Motor vehicles': 'realx', 'Nondurable goods': 'realc'}

model_names = dict()
model_names['Real'] = ['exp', 'realx', 'realc', 'pdur', 'r']
model_names['Nominal'] = ['nomexp', 'nomx', 'nomc', 'pdur', 'r']

Tmax = 12
modeldates = pd.date_range(start='2008-01-01', freq='MS', periods=Tmax)

ratelist = ['r','pi','rr','rk']

# graph defaults
plt.rc('font', size=12)
new_prop_cycle = (cycler('color', ['k','tab:blue','tab:purple','tab:green','r']) + 
                  cycler('linestyle', ['-','-','-','-','--']))
new_prop_cycle_psmj = (cycler('color', ['k','tab:purple','tab:green','r']) + 
                       cycler('linestyle', ['-','-','-','--']))
new_prop_cycle_model = (cycler('color', ['tab:blue','tab:purple','tab:green','r']) + 
                        cycler('linestyle', ['-','-','-','--']))
new_prop_cycle_orw = (cycler('color', ['k','tab:blue','r']) + 
                      cycler('linestyle', ['-','-','--']))                           
plt.rc('axes', prop_cycle=new_prop_cycle)

# ------------------------------------------------------------------------
# Data for counterfactual
# ------------------------------------------------------------------------

dfreal = pd.read_parquet('../input/pceexpenditurereal.parquet')
dfnom = pd.read_parquet('../input/pceexpenditure.parquet')
dfprice = pd.read_parquet('../input/pceexpenditureprice.parquet')

# relative price
dfrelprice = pd.read_stata('../input/motorvehicleprice.dta').rename(columns={'date': 'TimePeriod', 'newvehiclesresearchrelp': 'Relative Durable Price'}).set_index('TimePeriod')

# nominal interest rate
dfffr = pd.read_stata('../input/freddata_for_forecasting.dta').rename(columns={'mdate': 'TimePeriod', 'ffr': 'Federal Funds Rate'}).set_index('TimePeriod').loc[:, 'Federal Funds Rate']

dforg = dict()
for df, defl in zip([dfreal, dfnom],['Real','Nominal']):
    vars_to_plot = ['Personal consumption expenditures', 'Motor vehicles and parts', 'Nondurable goods', 'Services'] #,
    df = df['DataValue']
    df = df.loc[vars_to_plot].swaplevel().unstack()
    # combine nondurables and services
    df['Nondurable goods and services'] = df['Nondurable goods'] + df['Services']
    vars_to_plot.remove('Services')
    vars_to_plot.remove('Nondurable goods')
    vars_to_plot.append('Nondurable goods and services')
    # add relative durable price
    df = df.merge(dfrelprice, left_index=True, right_index=True)
    df = df.merge(dfffr, left_index=True, right_index=True)
    vars_to_plot.append('Relative Durable Price')
    vars_to_plot.append('Federal Funds Rate')
    df = df[vars_to_plot]
    dforg[defl] = df.copy()

# normalize to 2008 m1 nominal spending
for var in dforg['Real'].columns:
    dforg['Real'][var] = dforg['Real'][var] / dforg['Real'].loc['2008-01-01',var] * dforg['Nominal'].loc['2008-01-01',var]


# forecast
dffc = pd.read_stata('../input/forecasts.dta').rename(columns={'mdate': 'TimePeriod'}).set_index('TimePeriod')
dffcCI = pd.read_stata('../input/forecasts_lowerCI.dta').rename(columns={'mdate': 'TimePeriod'}).set_index('TimePeriod')

dffcset = {'Real': dffc['lrconsforA'].rename('Pessimistic Forecast'), 'Nominal': (dffc['lrconsforA']*dffc['lpconsforA']).rename('Pessimistic Forecast')}    
for defl in ['Real','Nominal']:
    dffcset[defl] =  dffcset[defl] / dffcset[defl].loc['2008-01-01'] * dforg[defl].loc['2008-01-01', 'Personal consumption expenditures']
dffcsetCI = {'Real': dffcCI['lrconslc2forA'].rename('Pessimistic Forecast: Lower CI'), 'Nominal': (dffcCI['lrconslc2forA']*dffcCI['lpconslc2forA']).rename('Pessimistic Forecast: Lower CI')}    
for defl in ['Real','Nominal']:
    dffcsetCI[defl] =  dffcsetCI[defl] / dffcsetCI[defl].loc['2008-01-01'] * dforg[defl].loc['2008-01-01', 'Personal consumption expenditures']
#%%
# add names (manual)
episodes = {'Date': ['Jan-Apr 2020', 'Jan-Apr 1980', 'Aug-Nov 1974', 'Apr-Jul 1960', 'Jul-Oct 2005', 'Sep-Nov 2008'],
            'Episode': ['COVID lockdowns', 'Credit controls, Volcker', 'prior spike up', 'prior spike up', 'prior spike up', 'Lehman Collapse']}
dates = pd.to_datetime(['2020-04-01','1980-04-01','1974-11-01','1960-07-01','2005-10-01','2008-11-01'])
dfnames = pd.DataFrame(data = episodes, index=dates)

dfdecline = dict()
for var, cutoff in zip(vars_to_plot[:2],[-1.4,-24]):
    # calculate growth rates
    pcenom = dfnom.loc[var,'DataValue']
    pceprice = dfprice.loc[var,'DataValue']

    pcereal = pcenom / pceprice

    # 3-month growth rates
    pcegrowth = pcereal.pct_change(periods = 3) * 100

    # largest 3-month declines
    largestdeclines = pcegrowth[pcegrowth<=cutoff]
    largestdeclines.loc['2008-11-01'] = pcegrowth.loc['2008-11-01']
    
    # strip out adjecent dates
    dates = largestdeclines.index
    maxdeclines = []
    for date in dates:
        # find adjecent dates
        currdates = dates[abs(dates - date).days < 300]
        currdecline = pcegrowth[currdates].idxmin()

        if currdecline not in maxdeclines:
            maxdeclines.append(currdecline)

    # table of max declines
    dfdeclinetemp = np.round(-pcegrowth[maxdeclines].sort_values(ascending=True).rename('Decline'),1)

    dfdeclinetemp = dfdeclinetemp.to_frame().merge(dfnames, how='left', left_index=True, right_index=True)

    dfdecline[var] = dfdeclinetemp[['Date', 'Episode', 'Decline']]

print(dfdecline)
#%%

# ------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------

def mpctable(data, table_dict, mpcset, simultitle, savename):

    GEset = data[list(mpcset)[0]]['MPCs'].keys()

    f = open('../output/' + savename + '.tex', 'w')
    f.write('\\begin{table}[htbp]\n')
    f.write('\\caption{General Equilibrium Marginal Propensity to Consume: ' + simultitle + '}\n') 
    f.write('\\begin{tabularx}{\\columnwidth}{@{\\hskip\\tabcolsep\\extracolsep\\fill}*{' + str(int(len(table_dict)*len(GEset))) + '}{S[table-format=1.2]}}\\toprule\n') 
    for i, title in enumerate(table_dict.keys()):
        title = '\\multicolumn{2}{c}{' + title + '}'
        if 1+i<len(table_dict):
            f.write(title + ' & ')
        else:
            f.write(title + ' \\\\  \n')
    for i, title in enumerate(table_dict.keys()):
        for j, GE in enumerate(GEset):
            title = '\\multicolumn{1}{c}{' + GE + '}'
            if 1+i==len(table_dict) and 1+j == len(GEset):
                f.write(title + ' \\\\ \\midrule \n')
            else:
                f.write(title + ' & ')
    for mpcval in mpcset:
        for i, var in enumerate(table_dict.values()):
            for j, GE in enumerate(GEset):
                f.write(f"{data[mpcval]['MPCs'][GE][var]['12 months cumulative']:.2f}")
                if 1+i==len(table_dict) and 1+j == len(GEset):
                    f.write('\\\\ \n')
                else:
                    f.write(' & ')
    f.write('\\bottomrule\\end{tabularx}\n') 
    f.write('\\label{tab:' + savename + '}\n')    
    f.write('\\end{table}\n')    
    f.close()


# ------------------------------------------------------------------------
# Tables and Figures
# ------------------------------------------------------------------------

#%%
for parasetname, parasettitle in zip(parasetnames, parasettitles):

    output = pickle.load( open( '../output/' + parasetname + '.pkl', 'rb' ) )
    parameter_set = output.keys()

    for parameter in parameter_set:

        mpcset = output[parameter].keys()
        info_set = output[parameter][list(mpcset)[0]]['shock'].keys()
        GEset = output[parameter][list(mpcset)[0]]['MPCs'].keys()

        # calculate MPCs:
        for info in info_set:
            for GE in GEset:
                for mpcval in mpcset:
                    for var in table_dict.values():
                        deltac = output[parameter][mpcval]['shock'][info][GE][var].sum() * output[parameter][mpcval]['steady'][var]
                        deltat = -output[parameter][mpcval]['shock'][info][GE]['t'].sum() * output[parameter][mpcval]['steady']['t']
                        output[parameter][mpcval]['MPCs'][GE][var]['12 months cumulative'] = deltac / deltat

        savetitle = 'mpcs' + parameter + parasetname
        mpctable(output[parameter], table_dict, mpcset, parasettitle, savetitle)
        
        if parasetname == 'baseline':
            savetitle = 'mpcs' + parameter + parasetname + 'PSMJ'
            mpctable(output[parameter], {'PCE': 'exp'}, list(mpcset)[1:], parasettitle, savetitle)

        elif parasetname == 'durableprice' or parasetname == 'lesselastic':
            savetitle = 'mpcs' + parameter + parasetname + 'ORW'
            mpctable(output[parameter], {'PCE': 'exp'}, {list(mpcset)[0]}, parasettitle, savetitle)            

    # create figures
    for parameter in parameter_set:
        for info in info_set:
            for GE in GEset:
                for defl in model_names.keys():
                    for varname, model_name, var in zip(plot_names, model_names[defl], vars_to_plot):
                        
                        if var=='Relative Durable Price' and defl=='Nominal':
                            continue

                        df = dforg[defl][var].copy().rename('Data').to_frame()

                        for mpcval in mpcset:
                            # if var=='Relative Durable Price' and parasetname=='durableprice':
                            #     stop
                            if model_name not in output[parameter][mpcval]['shock'][info][GE].keys():
                                output[parameter][mpcval]['shock'][info][GE][model_name] = np.zeros_like(output[parameter][mpcval]['shock'][info][GE]['exp'])
                            dfmodel = pd.DataFrame.from_dict(output[parameter][mpcval]['shock'][info][GE][model_name])

                            dfmodel.columns = [mpcval]
                            dfmodel = dfmodel.rename(index=dict(zip(range(Tmax), modeldates)))
                            df = df.merge(dfmodel, how='left', left_index=True, right_index=True)
                            
                            if model_name in ratelist:
                                df[mpcval] = df['Data']  - df[mpcval] * 100
                            else:
                                df[mpcval] = df['Data'] * ( 1 - df[mpcval] )

                        dfplot = df.loc['2008-01-01':'2008-12-01', :]

                        # if var=='Relative Durable Price':
                        #     stop

                        if varname=='PCE':
                            dfplotset = [dfplot, dfplot.merge(dffcset[defl], left_index=True, right_index=True),dfplot.merge(dffcsetCI[defl], left_index=True, right_index=True)]
                        #elif varname == 'PCE_CI':
                        #    dfplotset = [dfplot, dfplot.merge(dffcsetCI[defl], left_index=True, right_index=True)]
                        else:
                            dfplotset = [dfplot]

                        if var in vars_to_plot[:2] and GE=='GE':
                            dftab =  np.round(- (df.loc['2008-07-01', list(mpcset)] / df.loc['2008-04-01', list(mpcset)] - 1) * 100 , 1)

                            dftab = dftab.to_frame(name='Decline').reset_index().rename(columns={'index': 'Episode'})
                            dftab['Date'] = ''
                            
                            dfexp = pd.concat([dfdecline[var], dftab], axis=0, ignore_index=True).sort_values(by=['Decline'], ascending=False)
                            
                            dfexp.to_latex('../output/' + defl + '_' + varname.replace(' ','') + '_drop_' + GE + '_' + parasetname + '.tex', index=False)

                        # full plot
                        for dfplot, fc in zip(dfplotset, ['','fc','fclow']):
                            dfplotc = dfplot.copy()
                            # for dfpl, colors, psmj in zip((dfplotc, dfplotc.drop([list(mpcset)[0]], axis=1), dfplotc.drop(list(mpcset)[1:], axis=1)),(new_prop_cycle,new_prop_cycle_psmj,new_prop_cycle_orw),('','PSMJ','ORW')):
                            # for dfpl, colors, psmj in zip(dfplotc,new_prop_cycle,''):
                            dfpl = dfplotc
                            colors = new_prop_cycle
                            psmj = ''
                            fig, ax = plt.subplots()
                            # plt.rc('font', size=12) #controls default text size
                            ax.set_prop_cycle(colors)
                            # ax.set_color_cycle(['k','tab:blue','tab:orange','tab:green'])
                            # plt.plot(dfpl, linewidth=3)
                            for col in dfpl.columns:
                                if col=='Data':
                                    order = 2.5
                                else:
                                    order = 2
                                plt.plot(dfpl[col], linewidth=3, zorder=order)
                            plt.xlim(datetime.datetime(2008,1,1), datetime.datetime(2008,12,1))
                            # ax.xaxis.set_major_locator(mdates.MonthLocator())
                            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
                            plt.xlabel('2008')
                            # plt.xticks(rotation=45)
                            if varname == 'PCE' and defl == 'Real':
                                plt.ylim(790,840)
                            elif varname == 'PCE' and defl == 'Nominal':
                                plt.ylim(810,850)
                            elif varname == 'Motor vehicles' and defl == 'Real':
                                plt.ylim(0,35)
                            elif varname == 'Motor vehicles' and defl == 'Nominal':
                                plt.ylim(0,35)
                            if var=='Relative Durable Price':
                                plt.ylim(92,102)
                                plt.ylabel('Index, Jan. 2008 = 100')
                            elif var=='Federal Funds Rate':
                                plt.ylabel('Percent Annualized')
                            else:
                                plt.ylabel('Billion $')
                            # plt.title(defl + ' ' + varname)
                            plt.legend(list(dfpl.columns), loc='lower left')
                            plt.tight_layout(pad=0)
                            # plt.savefig('../output/mpc' + parameter + defl + varname.replace(' ','') + GE + info + parasetname + '.eps')
                            plt.savefig('../output/' + defl + '_' + varname.replace(' ','') + fc + '_' + GE + '_' + parasetname + psmj + '.eps')

                            # plt.title(defl + ' ' + varname)
                            ax.legend_ = None
                            # plt.savefig('../output/mpc' + parameter + defl + varname.replace(' ','') + GE + info + parasetname + '.eps')
                            plt.savefig('../output/' + defl + '_' + varname.replace(' ','') + fc + '_' + GE + '_' + parasetname + psmj + '_nolegend.eps')
                            # plt.show()
                            plt.close()
                        
                        # small plot
                        # dfplotsmall = dfplot[['Data'] + plotlist].copy()
                        # plt.plot(dfplotsmall)
                        # plt.xlim(datetime.datetime(2008,1,1), datetime.datetime(2008,12,1))
                        # plt.legend(list(dfplotsmall.columns), loc='lower left')
                        # plt.title(defl + ' ' + varname + ': ' + GE + ' ' + info)
                        # plt.savefig('../output/mpc' + parameter + defl + varname.replace(' ','') + GE + info + parasetname + 'small.eps')
                        # plt.show()
                        # plt.close()
                        
                # save all variables to excel:
                for mpcval in output[parameter].keys():
                    dfirf = pd.DataFrame.from_dict(output[parameter][mpcval]['shock'][info][GE])
                    dfirf = dfirf.rename(index=dict(zip(range(Tmax), modeldates)))


                    dfirf.to_excel('../output/irfs' + parameter + 'mpc' + mpcval[-2:] + GE + info + parasetname + '.xlsx')

                # plot durable price IRF
                if GE == 'GE':
                    for plotvar, plotname in zip(['pdur','pi','r'],['Relative Durable Price','Inflation','Nominal Interest Rate']):
                        df = pd.DataFrame(index = modeldates)
                        for mpcval in output[parameter].keys():
                            dfirf = pd.DataFrame.from_dict(output[parameter][mpcval]['shock'][info][GE][plotvar] * 100)
                            dfirf.columns = [mpcval]
                            dfirf = dfirf.rename(index=dict(zip(range(Tmax), modeldates)))
                            df = df.merge(dfirf, how='left', left_index=True, right_index=True)

                        fig, ax = plt.subplots()
                        ax.set_prop_cycle(new_prop_cycle_model)
                        plt.rc('font', size=12) #controls default text size
                        plt.plot(df)
                        plt.xlim(datetime.datetime(2008,1,1), datetime.datetime(2008,12,1))
                        # plt.xticks(rotation=45)
                        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
                        plt.xlabel('2008')
                        if plotvar not in ['r','pi','rr','rk']:
                            plt.ylabel('% Deviation from SS')
                            plt.legend(list(df.columns), loc='upper left')
                        elif plotvar in ['r']:
                            plt.ylabel('Basis Points')
                            plt.legend(list(df.columns), loc='upper left')
                        else:
                            plt.ylabel('Percentage Point Deviation from SS')
                            plt.legend(list(df.columns), loc='lower left')
                        
                        # plt.title(defl + ' ' + varname)
                        plt.tight_layout()
                        # plt.savefig('../output/mpc' + parameter + defl + varname.replace(' ','') + GE + info + parasetname + '.eps')
                        plt.savefig('../output/' + plotvar + '_' + GE + '_' + parasetname + '.eps')
                        # plt.show()
                        plt.close()
                    
                
                    



# ------------------------------------------------------------------------
# Output parameters
# ------------------------------------------------------------------------

        ss = output[parameter][list(mpcset)[0]]['steady']

        tablenames = {
            'beta': ['$\\beta$', 'Subjective discount factor'],
            'psi': ['$\\psi$', 'Weight on durable service flow'],
            'sigma': ['$\\sigma$', 'Utility curvature on nondurable consumption'],
            'sigmad': ['$\\sigma^d$', 'Utility curvature on durable service flow'],
            'thetad': ['$\\theta^d$', 'Calvo parameter on durable adjustment'],
            'eta': ['$\\eta$', 'Durable operating cost'],
            'nu': ['$\\nu$', 'Weight on disutility of labor'],
            'phi': ['$\\phi$', 'Inverse of the Frisch elasticity of labor supply'],
            'gamma': ['$\\gamma$', 'Fraction of Hand-to-Mouth consumers'],
            'mpx': ['$\\vartheta$', 'Hand-to-Mouth MPC on durables'],
            'deltad': ['$\\delta^d$', 'Depreciation of durable consumption goods'],
            'alpha': ['$\\alpha$', 'Exponent on private capital in production function'],
            'delta': ['$\\delta$', 'Depreciation of private capital'],
            'kappa': ['$\\kappa$', 'Investment adjustment cost parameter'],
            'delta1': ['$\\delta_1$', 'Parameter on linear term of capital utilization cost'],
            'delta2': ['$\\delta_2$', 'Parameter on quadratic term of capital utilization cost'],
            'zeta': ['$\\zeta$', 'Inverse relative supply elasticity of durable goods'],
            'mup_ss': ['$\\mu_p$', 'Steady-state price markup'],
            'muw_ss': [ '$\\mu_W$', 'Steady-state wage markup'],
            'theta': ['$\\theta_p$', 'Calvo parameter on price adjustment'],
            'thetaw': ['$\\theta_W$','Calvo parameter on wage adjustment'],
            'eps': ['$\\epsilon_p$', 'Elasticity of substitution between types of goods'],
            'epsw': ['$\\epsilon_W$', 'Elasticity of substitution between types of labor'],
            'gyfrac': ['$gy$', 'Steady-state share of total govt spending to GDP'],
            'phib': ['$\\phi_b$', 'Debt feedback coefficient in fiscal rule'],
            'rhor': ['$\\rho_{r}$', 'Monetary policy interest rate smoothing'],
            'phipi': ['$\\phi_{\\pi}$', 'Monetary policy response to inflation'],
            'phigap': ['$\\phi_{gap}$', 'Monetary policy response to the output gap'],
        }
        tabletitles = ['Parameter', 'Value', 'Description']
        tableformat = ['l','S[table-format=2.3]','l']

        # create smaller version
        presentation = ['sigma','phi','gamma','thetad','sigmad','mpx','psi','deltad','zeta','phib']
        tablenames_presentation = dict()
        for para in presentation:
            tablenames_presentation[para] = tablenames[para]

        # nondurable version
        nondur = ['beta','psi','sigma','nu','phi','gamma','alpha','kappa','delta1','delta2','mup_ss','muw_ss','theta','thetaw','eps','epsw','gyfrac','phib','rhor','phipi','phigap']
        tablenames_nondur = dict()
        for para in nondur:
            tablenames_nondur[para] = tablenames[para]

        # durable version
        dur = ['sigmad','thetad','deltad','eta','mpx','zeta']
        tablenames_dur = dict()
        for para in dur:
            tablenames_dur[para] = tablenames[para]

        for table, tablesave in zip([tablenames, tablenames_presentation, tablenames_nondur, tablenames_dur],['calibration','calibration_small','calibration_nondur','calibration_dur']):
            f = open('../output/' + tablesave + parasetname + parameter + '.tex', 'w')
            f.write('\\begin{table}[htbp]\n')
            if tablesave!='calibration_small':
                f.write('\\caption{Baseline Calibration of the Model}\n') 
            f.write('\\begin{tabular}{' + ''.join(tableformat) + '}\\toprule\n') 
            for i, title in enumerate(tabletitles):
                if tableformat[i][0]=='S':
                    title = '\\multicolumn{1}{l}{' + title + '}'
                if 1+i<len(tabletitles):
                    f.write(title + ' & ')
                else:
                    f.write(title + ' \\\\ \\midrule \n')
            for param, text in table.items():
                if ss[param]==0:
                    continue
                if param in ['gamma','thetad', 'mpx']:
                    f.write(text[0] + ' & \\text{varies}')
                else:
                    f.write(text[0] + ' & ' + str(round(ss[param], 3)) )
                for col in text[1:]:
                    f.write(' & ' + col)
                f.write('\\\\ \n')
            f.write('\\bottomrule\\end{tabular}\n') 
            f.write('\\begin{minipage}{\hsize} \\rule{0pt}{9pt} \\footnotesize \n')
            f.write('Notes: The model is calibrated at a monthly frequency. ')
            if 'gamma' in table.keys():
                f.write('The parameter $\\gamma$ is calibrated to either 0.3, 0.5, or 0.9, which corresponds to the aggregate MPC in the model. ')
            if 'thetad' in table.keys():
                f.write('The parameter $\\theta^d$ is calibrated such that for each value of $\\gamma$ to model replicates our empirical targets for the short-term interest elasticity of durable demand. ')
                f.write('For example, when $\\gamma=' + str(round(ss['gamma'], 2)) + '$, then $\\theta^d=' + str(round(ss['thetad'], 3)) + '$. ')
            if 'mpx' in table.keys():
                f.write('The parameter $\\vartheta$ is calibrated  to match an overall MPC on motor vehicles of 0.3 when $\\gamma=0.3$ and of 0.4 when $\\gamma=0.5$ or $\gamma=0.9$. ')
                f.write('This yields  $\\vartheta=' + str(round(0.3/0.3, 2)) + '$ when $\\gamma=0.3$,  $\\vartheta=' + str(round(0.4/0.5, 2)) + '$ when $\\gamma=0.5$, and  $\\vartheta=' + str(round(0.4/0.9, 2)) + '$ when $\\gamma=0.9$. ')
            f.write('See the text for details. ')
            f.write('\\end{minipage}\n')
            f.write('\\label{tab:' + tablesave + '}\n')    
            f.write('\\end{table}\n')    
            f.close()

# %%
