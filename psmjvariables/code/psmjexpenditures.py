#%%
import pandas as pd
import numpy as np

                           
# expenditure variables to be used in regression
# for the first set we will also add the log of expenditure
keep_variables_log = ['FOODBEVS','SNDEXP','NDEXP','TOTEXP2','TOTEXP_NMV'] 
keep_variables_nolog = ['CARTKN','CARTKU','NEWTRUCKS','NEWPCARS','USED_CAR_PCE','DEXP','SEXP','EDUCA', 'MAINRP', 'VEHINS', 'ENTERT',  'VRNTLO',  'VEHFIN',  'SHELT', 'FINCHRG', 'MORTGAGEINT']
keep_variables_all = keep_variables_log + keep_variables_nolog

# read in monthly expenditure files
dfmtbi = pd.read_parquet('../input/mtbimonthly.parquet', 
                         columns=keep_variables_all)

#adds all psmj categories
dfmtbi_add = pd.read_parquet('../input/mtbimonthly.parquet')
dfmtbi_add = dfmtbi_add.add_prefix('psmj_')

# add interview date to dataframe and set index to household - interview date
dffmli = pd.read_parquet('../input/fmliquarterly.parquet', columns = ['NEWID','CUID','INTDATE'])




print('The Number of observations in the fmli file is:')
print(str(dffmli.shape[0])) # There are 82909 in 2007-2009 sample and 82897 in the trimmed 1996-2019 sample

print('The Number of observations in the mtbi file is:')
print(str(dfmtbi.shape[0])) # There are 234789 in 2007-2009 sample and 249188 in the trimmed 1996-2019 sample


#%%
#Merge Expenditure Data with family characteristics
dfmtbi = dfmtbi.merge(dffmli, how='left', right_index=True, left_index=True)
dfmtbi = dfmtbi.merge(dfmtbi_add, how='left', right_index=True, left_index=True)


# ------------------------------------------------------------------------
# Loop over aggregation
# - Creates logs and differences for monthly or 3-month expenditures
# ------------------------------------------------------------------------

# frequency at which data is recorded  
timefreq = 'MS'

# could run two versions, monthly or interview level
for expenditure_frequency in ['monthly', 'interview']: 

    # existing data already at monthly frequency
    if expenditure_frequency=='monthly':
        timevar = 'DATE'
        timeshift = 1           # expenditures are one month apart
        
        # set index to household ID and expenditure month
        indexset = ['CUID', 'NEWID', timevar, 'INTDATE']
        df = dfmtbi.reset_index().set_index(indexset) 


    # for interview level we aggregate and merge interview date
    elif expenditure_frequency=='interview':   
        timevar = 'INTDATE'
        timeshift = 3           # interviews are three months apart

        # aggregate to interview level 
        indexset = ['CUID', 'NEWID', timevar]
        df = dfmtbi.groupby(indexset).sum()


    # create log of each expenditure variable
    # log will create -inf when arguemnt is 0, so reaplce with nan
    dflog = df[keep_variables_log].apply(np.log).replace((np.inf, -np.inf), np.nan)

    df = df.merge(dflog.add_prefix('l_'), how='inner', left_index=True, right_index=True)

    # next set of lines first difference the data
    vars_to_difference = df.columns 

    # have date as sole index so we can use differences
    dfshift = df[vars_to_difference].reset_index().set_index(timevar)

    # shift columns by timeshift periods of frequency timefreq
    dfshift = dfshift.groupby('CUID').shift(periods=timeshift, freq=timefreq)
    dfshift = dfshift.reset_index().set_index(indexset)

    # substract shifted data to create difference (must be on same index again)
    dfdiff = df[vars_to_difference] - dfshift[vars_to_difference]

    # rename all the columns
    dfdiff = dfdiff.add_prefix('d_')
    dflag = dfshift.add_prefix('LAG_')

    # # merge with the original dataset keeping only the row observations we had originally
    df = df.merge(dfdiff, how='left', left_index=True, right_index=True)
    df = df.merge(dflag, how='left', left_index=True, right_index=True)
    
    df.to_parquet('../output/psmjexpenditures' + expenditure_frequency + '.parquet')
    print(str(df.shape[0])) # There are 234789 in 2007-2009 sample and 249188 in the trimmed 1996-2019 sample

#%%
# ------------------------------------------------------------------------
# Testing aggregation
# ------------------------------------------------------------------------
if __name__=='__main__':
    df = pd.read_parquet('../output/psmjexpendituresmonthly.parquet', columns = keep_variables_all)
    df = df.groupby(['CUID', 'NEWID']).sum()    

    dfint = pd.read_parquet('../output/psmjexpendituresinterview.parquet', columns = keep_variables_all)
    dfint = dfint.add_suffix('_int')
    # df = dfint.groupby(['CUID', 'NEWID']).sum()

    df = df.merge(dfint, how='left', left_index=True, right_index=True)

    for var in keep_variables_all:
        print('Aggregation Residuals for ' + var + ':')
        print((df[var]-df[var + '_int']).abs().sum())
        print((df[var].isna()^df[var + '_int'].isna()).abs().sum())


# %%
