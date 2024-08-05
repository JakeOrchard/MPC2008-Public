#%%
import pandas as pd
import numpy as np
from cohortweights import cohortweights
from staggered_treatment_MC import staggered_treatment_MC

# TODOHERE:
# Add control variables?

#%%
# ------------------------------------------------------------------------
# Monte Carlo Parameters
# ------------------------------------------------------------------------

np.random.seed(2021) 

S = 10 #100
T = 4
TC = 2
Nset = [2,8,20]
treatvals = [-1,1]
treattime = [0,1]

# print((np.array(Nset[0:len(treattime)]) * np.array(treatvals)).sum() / np.array(Nset[0:len(treattime)]).sum())


# ------------------------------------------------------------------------
# Simulate dataset and Run Regressions
# ------------------------------------------------------------------------


weight_est = []

for s in range(S):

    print(s)

    df = staggered_treatment_MC(T=T, TC=TC, Nset=Nset, treatvals=treatvals, treattime=treattime)
    
    # set index to household to facilitate household level calculations
    df.set_index('HHind', inplace=True)

    # First date by HH
    df = df.merge(df['Date'].groupby(['HHind']).min().rename('FirstDate'), left_index=True, right_index=True, how='left')

    # First treatment by HH
    df = df.merge(df.loc[df['D']==1, 'Date'].groupby(['HHind']).min().rename('FirstTreat'), left_index=True, right_index=True, how='left')
    
    # If treatment is missing replace with Date + Large Number
    nevertreated = pd.isna(df['FirstTreat'])
    df.loc[nevertreated, 'FirstTreat'] = df.loc[nevertreated, 'FirstDate'] + 10 ** 6
    
    # Cohort is group by first date and first treatment
    df['Cohort2'] = df.groupby(['FirstDate','FirstTreat']).ngroup()

    # construct relative time indicator
    df['RelTime2'] = df['Date'] - df['FirstTreat']

    # Then set panel dimension
    df = df.reset_index().set_index(['HHind','Date'])


    weight_est0, dfweights = cohortweights(df, 
                                           outcome='Y', 
                                           treatment='D', 
                                           relativetime='RelTime2', 
                                           cohort='Cohort2')

    weight_est.append(weight_est0)

print(weight_est)
# %%
