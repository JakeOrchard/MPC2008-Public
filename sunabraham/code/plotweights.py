#%%
import pandas as pd
import numpy as np
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import os, sys
rootpath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(rootpath)
from utilities.mpc_utilities import savetexscalar

# TODO:
# Treatments should be same, also very negative later in the sample does it make sense?
# sould i look at different relative times for leads and lags?\

for file in sys.argv[1:]:
# for file in [../output/weightsdtotexp2lags0insampleinterview.parquet]: #

    print(file)

    # ------------------------------------------------------------------------
    # Load Data
    # ------------------------------------------------------------------------

    # filename = 'weights' + outcome + 'lags' + str(lags) + sample + frequency 
    # filename = filename.replace('_','').replace(' ','').lower()

    dfweights = pd.read_parquet(file)  

    filename, extension = os.path.splitext(os.path.basename(file))

    title = 'get title'

    # extract date variable
    for var in dfweights.index.names:
        if var.endswith('DATE'):
            datevar = var

    dfsum = dfweights['Contribution'].groupby('Variable').sum()                    

    for var in dfsum.index:
        if np.abs(dfsum[var])>10:
            decimals = 1
        else:
            decimals = 2

        name= 'scalar' + filename + var.replace(' ','').lower()
        validname = re.sub(r'[0-9]|\_+', '', name)        

        savetexscalar(value = dfsum[var], 
                    name= validname, 
                    decimals = decimals)
    
    # ------------------------------------------------------------------------
    # Aggregate to treatment and control groups
    # ------------------------------------------------------------------------

    bytreat = dfweights
    bytreat['Treat Status'] = np.sign(bytreat.index.get_level_values(level='RELATIVE TIME'))
    bytreat = bytreat.groupby([datevar, 'Variable', 'Treat Status']).sum()
    bytreat['Treatment'] = bytreat['Contribution'] / bytreat['Weights']
    bytreat = bytreat.drop('Cohort Treatment', axis=1)

    bytreat = bytreat.loc['2008-06-01':]
    bytreat = bytreat.unstack().fillna(0).stack()

    # print(bytreat)

    # Treatments should be same, also very negative later in the sample does it make sense?

    # ------------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------------

    treatment_list = bytreat.index.get_level_values(level='Variable').unique()
    date_list = bytreat.index.get_level_values(level=datevar).unique()

    width = 8
    delta = pd.Timedelta(width, unit='D')
    
    for treat in treatment_list:

        for plotvar in ['Treatment','Contribution','Weights']:
            fig, ax = plt.subplots()

            for group, timevalue in {'Control': -1, 'Treatment': 0, 'PastTreat': 1}.items():
                ax.bar(date_list + timevalue * delta, bytreat.loc[pd.IndexSlice[:, treat, timevalue],plotvar] , width, label=group)

            # rects1 = ax.bar(date_list - delta, bytreat.loc[pd.IndexSlice[:, var,-1],plotvar] , width, label='Control')
            # rects2 = ax.bar(date_list, bytreat.loc[pd.IndexSlice[:, var, 0],plotvar], width, label='Treatment')
            # rects3 = ax.bar(date_list + delta, bytreat.loc[pd.IndexSlice[:, var, 1],plotvar] , width, label='PastTreat')

            ax.legend()

            # Add some text for labels, title and custom x-axis tick labels, etc.
            if plotvar!='Weights':
                ax.set_ylabel('\$')
            ax.set_xticks(date_list)
            
            # Minor ticks every month.
            fmt_month = mdates.MonthLocator()
            ax.xaxis.set_major_locator(fmt_month)
            fig.autofmt_xdate()

            ax.set_title(plotvar + ': ' + title)

            figname = filename.replace('weights','')
            suffix = (treat + plotvar).replace(' ','').lower()

            plt.savefig('../output/' + figname + suffix + '.png')
            
            # plt.close()


    # overall = bytreat.groupby(datevar).agg({'Contribution': 'sum', 'Weights': 'max'})
    # overall['Treatment'] = overall['Contribution'] / overall['Weights']    

    # for plotvar in ['Treatment','Contribution']:
    #     fig, ax = plt.subplots()
    #     rects1 = ax.bar(x, overall.loc[:,plotvar] , width, label='Treatment Effect')
    #     ax.legend()  

    #     # Add some text for labels, title and custom x-axis tick labels, etc.
    #     if plotvar!='Weights':
    #         ax.set_ylabel('\$')
    #     ax.set_xticks(x)
        
    #     # Major ticks every month.
    #     fmt_month = mdates.MonthLocator()
    #     ax.xaxis.set_major_locator(fmt_month)
    #     fig.autofmt_xdate() 

    #     ax.set_title(plotvar + ': ' + title)

    #     plt.savefig('../output/' + outcome + sample + plotvar + 'overall.png')

    #     plt.close()



# # %%
# # cannot cumulate with 3 month averaged expenditure data
# outcomelvl = outcome.replace('d_','')

# dfagg = df[[outcomelvl,'FINLWT21']]
# dfagg[outcomelvl] = dfagg[outcomelvl]*dfagg['FINLWT21']
# dfagg = dfagg[[outcomelvl,'FINLWT21']].groupby('INTDATE').sum()
# dfagg[outcomelvl] = dfagg[outcomelvl]/dfagg['FINLWT21']

# # cumulativenotreatment = bytreat.loc[pd.IndexSlice[:,False],'RT Treatment'].cumsum()
# # cumulativenotreatment = cumulativenotreatment.reset_index().set_index('INTDATE')

# # dfagg[outcomelvl + ' NO TREATMENT'] = dfagg[outcomelvl] + cumulativenotreatment['RT Treatment']

# print(dfagg)

# %%


    # outcomesets = {'d_CARTKN': 'Change in New Car Expenditure',
    #                     'd_TOTEXP2': 'Change in Total Expenditure',
    #                     'd_NDEXP': 'Change in Nondurable Expenditure'}