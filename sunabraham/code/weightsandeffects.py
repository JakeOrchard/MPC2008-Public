#%%
import pandas as pd
import numpy as np
from cohortweights import cohortweightsmult

# import os, sys
# rootpath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.append(rootpath)
# from utilities.mpc_utilities import savetexscalar

frequency_set = {'interview'} # , 'monthly'   monthly not working with lags right now
sample_set = {'INSAMPLE', 'INSAMPLE RBT'}
outcome_set = {'d_CARTKN', 'd_TOTEXP2', 'd_NDEXP'}
lag_set = [0, 1]
weights ='FINLWT21'  # None 
time_effects = True

# ------------------------------------------------------------------------
# Loop over interview and monthly sample
# ------------------------------------------------------------------------

for frequency in frequency_set:
    
    df = pd.read_parquet('../input/psmjsample' + frequency + '.parquet')

    # extract date variable
    for var in df.index.names:
        if var.endswith('DATE'):
            datevar = var

    # ------------------------------------------------------------------------
    # Create Cohort and Relative Time indicator
    # ------------------------------------------------------------------------
    # Cohort groups: when they enter the sample and when they are first treated 
    df['COHORT'] = df.groupby(['FIRSTDATE','FIRST RBT ' + datevar]).ngroup()  

    # construct relative time indicator in months
    df['RELATIVE TIME'] = ( df.index.get_level_values(level=datevar).year * 12 
                           - df['LAST RBT ' + datevar].dt.year * 12
                           + df.index.get_level_values(level=datevar).month
                           - df['LAST RBT ' + datevar].dt.month ) 

    cohort_var = ['COHORT']
    relativetime_var = ['RELATIVE TIME']
                               
    # for two-way FE regressions can have at most two levels
    df = df.droplevel('NEWID') 
    if frequency=='monthly':
        df = df.droplevel('INTDATE')        
        
    
    # ------------------------------------------------------------------------
    # Regression
    # ------------------------------------------------------------------------

    for lags in lag_set:


        # lag structure
        treat_vars = ['RBT INDICATOR']
        for lag in range(1,1+lags):
            treat_vars.append('LAG' + str(lag) + treat_vars[0])

        for outcome in outcome_set:
            
            if outcome.startswith('d_'):
                control_vars = ['AGE', 'd_PERSLT18', 'd_NUM_ADULTS']
                entity_effects = False  
            else:
                control_vars = ['AGE', 'PERSLT18', 'NUM_ADULTS']
                entity_effects = False  

        
            for sample in sample_set:

                print(frequency, lags, outcome, sample)

                # dfstata = df
                # # converts all True/False booleans to numeric
                # for col in dfstata.select_dtypes(include='bool').columns:
                #     dfstata[col] = dfstata[col].astype(int)

                # for col in dfstata.select_dtypes(include='object').columns:
                #     dfstata[col] = dfstata[col].multiply(1)
                #     dfstata[col] = dfstata[col].astype('float64')
                # dfstata.to_stata('../output/tempfile.dta')

                dfweights = cohortweightsmult(
                                            df.loc[df[sample]==1,], 
                                            outcome=[outcome], 
                                            treatment=treat_vars, 
                                            controls=control_vars,
                                            relativetime=relativetime_var, 
                                            cohort=cohort_var,
                                            time_effects=time_effects,
                                            entity_effects=entity_effects,
                                            weights=weights,
                                            )
                
                filename = 'weights' + outcome + 'lags' + str(lags) + sample + frequency
                filename = filename.replace('_','').replace(' ','').lower()

                dfweights.to_parquet('../output/' + filename + '.parquet')   

                # print(dfweights)    
                # bytreat = dfweights
                # bytreat['Treat Status'] = (bytreat.index.get_level_values(level='RELATIVE TIME')==0)
                # bytreat = bytreat.groupby([datevar, 'Variable', 'Treat Status']).sum()
                # bytreat['Treatment'] = bytreat['Contribution'] / bytreat['Weights']
                # bytreat = bytreat.drop('Cohort Treatment', axis=1)

                # print(bytreat)                                 

# %%
