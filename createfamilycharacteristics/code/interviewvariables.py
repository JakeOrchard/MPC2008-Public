#%%
import pandas as pd
import numpy as np

# ------------------------------------------------------------------------
# Loading and aggregating FMLI files
# ------------------------------------------------------------------------

fmli = pd.read_parquet('../input/fmliquarterly.parquet')


print('The Number of observations in the fmli file is (interviewvariables.py):')
print(str(fmli.shape[0])) # There are 82909 in 2007-2009 sample and 82897 in the trimmed 1996-2019 sample

#%%
# ------------------------------------------------------------------------
# Create interview variables from existing data
# ------------------------------------------------------------------------

# reset interview date
fmli['INTDATE'] = pd.to_datetime(fmli['INTDATE'])

# number of adults
fmli['NUM_ADULTS'] = fmli['FAM_SIZE'] - fmli['PERSLT18']

# age variable
# fix missing values in age 
# MOVE THIS TO PREPROCSSING STEP?
fmli['AGE2'] = fmli['AGE2'].replace(".","").apply(pd.to_numeric)

# average age if married, otherwise using household head age
fmli['AGE'] = fmli['AGE_REF']  

# for married use average age
fmli.loc[fmli['MARITAL1']==1, 'AGE']  = (fmli['AGE_REF'] + fmli['AGE2'] )/2



#Create Deciles of Income by interview date or pooled sample

fmli['INCD INTERVIEW'] = fmli.groupby(['INTDATE'])['FINCBTAX'].transform(
                     lambda x: pd.qcut(x, 10, labels = False,duplicates = 'drop'))

fmli['INCD POOL'] = pd.qcut(fmli['FINCBTAX'],10, labels = False,duplicates = 'drop')


if fmli['FINCBTAX'].isnull().sum() > 0:
    minmiss_inc = fmli[fmli['FINCBTAX'].isnull()].reset_index()['INTDATE'].min().strftime('%m/%d/%Y')
    maxmiss_inc = fmli[fmli['FINCBTAX'].isnull()].reset_index()['INTDATE'].max().strftime('%m/%d/%Y')
    print('ALERT: SOME HOUSEHOLDS BETWEEN ' + minmiss_inc + ' and ' + maxmiss_inc + ' MISSING INCOME INFORMATION')
    #CEX had another name for pre-tax income during missing dates. These dates are outside of
    #both the 2001 and 2008 rebate windows

#%%


#Log Income

fmli['LINCOME'] = fmli['FINCBTAX'].apply(np.log).replace((np.inf, -np.inf), np.nan)

# set new index: CUID and interview date
fmli = fmli.set_index(['CUID','INTDATE'])
#%%
# ------------------------------------------------------------------------
# First difference the data (if numeric)
# ------------------------------------------------------------------------
#%%
# interviews are three months apart
timeshift = 3
timefreq = 'MS'

# numeric columns
fmli_numeric = fmli.select_dtypes(include=np.number)

# move intdate as first level so we can unstack household level
fmli_numeric = fmli_numeric.reset_index().set_index('INTDATE')

# shift numeric columns (commented code is much slower)
# fmli_shifted = fmli_numeric.unstack().shift(periods=timeshift, freq=timefreq).stack()
fmli_shifted = fmli_numeric.groupby('CUID').shift(periods=timeshift, freq=timefreq)

# substract shifted data to create difference (must be on same index again)
fmli_numeric = fmli_numeric.reset_index().set_index(['CUID','INTDATE'])
fmli_diff = fmli_numeric - fmli_shifted

# rename all the columns
fmli_diff = fmli_diff.add_prefix('d_')
#%%
# merge with the original dataset keeping only the row observations we had originally
fmli = fmli.merge(fmli_diff, how='left', left_index=True, right_index=True)
print('here')
# ------------------------------------------------------------------------
# Save Output
# ------------------------------------------------------------------------

# save data
fmli.to_parquet('../output/interviewvariables.parquet')


#%%