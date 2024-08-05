# This file loads the CEX rebate module and processes it for merging with 
# the consumption expenditure files
#%%
import pandas as pd
import numpy as np
import os, sys
rootpath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(rootpath)
from utilities.mpc_utilities import aggregate_df, fill_missing_by_type


# ------------------------------------------------------------------------
# LOAD THE REBATE FILE
# ------------------------------------------------------------------------

df = pd.read_parquet('../input/rbt.parquet')

# get all households interview dates
intdate = pd.read_parquet('../input/fmliquarterly.parquet', columns = ['NEWID', 'CUID','RBTINTVIEW','INTDATE','INTNUM','FINLWT21'])

# get all households expenditure dates
hhdates = pd.read_parquet('../input/mtbimonthly.parquet', columns = ['DATE', 'NEWID'])

# merge the two and set index to household, date, interview date
hhdates = hhdates.merge(intdate, how='left', left_index=True, right_index=True)

hhdates = hhdates.reset_index().set_index(['CUID','DATE', 'NEWID','INTDATE'])

hhdates = hhdates.drop(['FINLWT21'], axis=1)
#%%
# ------------------------------------------------------------------------
# Create DATE and Rebate Number
# ------------------------------------------------------------------------

# frequency at which data is recorded  
timefreq = 'MS'

# number of rebates in a row
df['REBATES'] = 1

# Construct rebate year: has to besuch that the DATE occurs at or before to 
# interview date (multiple interviews, can cause rebate to appear to early)
df['RBTYR'] = np.where(df['RBTMO']>3,2008,2009) 

# DATE
df['DATE'] = pd.to_datetime(dict(year = df['RBTYR'], 
                                 month = df['RBTMO'], 
                                 day = 1))

#  Checks that no rebate received in April 2009
if df.loc[df['RBTMO']>=4, 'DATE'].max().year > 2008:
    print('REBATE RECEIVED AFTER MARCH 2009. INCONSISTENT WITH DATA DESCRIPTION')


# ------------------------------------------------------------------------
# Create Dummies from Categories
# ------------------------------------------------------------------------

df['RBT INDICATOR'] = True

# check or EFT indicators
df['CHECK'] = df['CHCKEFT']==1
df['EFT']   = df['CHCKEFT']==2

#Missing information on Check or EFT
if df['CHCKEFT'].isnull().sum() >0:
    print('REBATE SAMPLE MISSING DIST. INFO for ' +str(df['CHCKEFT'].isnull().sum()) + ' Rebate Observations' )

#%%
# Amt by check or EFT
df['RBTAMTCHECK'] = df['RBTAMT']*df['CHECK']
df['RBTAMTEFT'] = df['RBTAMT']*df['EFT'] 
     

# spend save indicator
df['SPEND'] = df['HOWUSED']==1
df['SAVE']  = (df['HOWUSED']==2) | (df['HOWUSED']==3)

# using qualitative answers
df['RBT SPEND'] = df['RBTAMT'] * df['SPEND']

# Drop variables we no longer use  
month_and_year = ['RBTMO','RBTYR','USDINTMO','USDINTYR']  #month and year variable
flags = ['RBTMO_','RBTAMT_','CHCKEFT_','HOWUSED_']   # data is either there or nan, flags do not add more info
categories = ['HOWUSED','CHCKEFT']

df = df.drop(month_and_year + categories + flags, axis=1)

# ------------------------------------------------------------------------
# Collapse at household - DATE level
# ------------------------------------------------------------------------
#%%
sum_vars = {'RBTAMT', 'RBT SPEND', 'REBATES', 'RBTAMTCHECK', 'RBTAMTEFT'}
max_vars = {'CHECK', 'EFT', 'RBT INDICATOR'}
min_vars = {'RBT DATE', 'RBT INTDATE'}

df = aggregate_df(df, 
                        agg_by={'NEWID', 'CUID', 'DATE'}, 
                        sum=sum_vars, 
                        max=max_vars)

# ------------------------------------------------------------------------
# Match DATE to interview date
# ------------------------------------------------------------------------
#%%
# merge interview date and number
df = df.merge(intdate[['INTDATE','INTNUM','FINLWT21']], how='left', left_index=True, right_index=True)

# report total sum of rebate
df['RBTWGHT'] = df['RBTAMT'] * df['FINLWT21']
print('Total Rebate Amount: ' + str(df['RBTWGHT'].sum()))

# drop weight ont total sum column
df = df.drop(['FINLWT21','RBTWGHT'], axis=1)

# we adjust the interview date next, so can no longer use NEWID as identifying index
# since it is tied to an interview date
df = df.reset_index().set_index(['CUID', 'DATE'])

# Distance between rebate (months)
df['DISTANCE RBT ITW'] = (
                            (df['INTDATE'].dt.year - df.index.get_level_values(level='DATE').year)*12 
                         +   df['INTDATE'].dt.month - df.index.get_level_values(level='DATE').month
                         )

# Check if rebate was received in the current reference period
df['WITHIN REF PERIOD'] = (df['DISTANCE RBT ITW']>0) & (df['DISTANCE RBT ITW']<=3)

# Baseline DATEs
df['MATCHED INTERVIEW'] = df.loc[df['WITHIN REF PERIOD']==True, 'INTDATE']

# Adjustment (1):
# For these households the rebate was received more than 3 months ago
# the appropriately matched interview is 3 months ago
# but the interview may not exist if the current interview is the first one
match_to_prev_intview = (df['DISTANCE RBT ITW'] == 4) & (df['INTNUM'] > 2)

df.loc[match_to_prev_intview, 'MATCHED INTERVIEW'] = df['INTDATE'] - pd.DateOffset(months=3) 

# Adjustment (2):
# check if DATE = interview date
# if so, set interview date is three month after DATE
# but first check if such an interview date exists
match_to_next_intview = (df['DISTANCE RBT ITW'] == 0) & (df['INTNUM'] < 5)

df.loc[match_to_next_intview, 'MATCHED INTERVIEW'] = df['INTDATE'] + pd.DateOffset(months=3) 
# df.loc[(df['DISTANCE RBT ITW'] == 0) & (df['INTNUM'] == 5), 'MATCHED INTERVIEW'] = pd.NaT 

# check that we matched all interviews unless we cannot match them to the next or 
# previous interview
assert ((df['MATCHED INTERVIEW'].isna()==False) | (match_to_next_intview.eq(0)) | (match_to_next_intview.eq(0))).all()

# keep only matched interviews:
df = df[df['MATCHED INTERVIEW'].isna()==False]


# rename interview since there is now a distinction when the rebate is reported and what interview we match to
df = df.drop(['INTDATE', 'INTNUM', 'NEWID', 'WITHIN REF PERIOD', 'DISTANCE RBT ITW'], axis=1)
df = df.rename({'MATCHED INTERVIEW': 'INTDATE'}, axis=1, errors="raise")


# These variables store the household rebate data
df = df.reset_index()
df['RBT DATE'] = df['DATE']
df['RBT INTDATE'] = df['INTDATE']

# put matched interview date in index for merging
df = df.set_index(['CUID', 'DATE', 'INTDATE'])

# ------------------------------------------------------------------------
# merge rebate data with all household-date observation and fill in missing values
# do this at household level, monthly, and interview frequency
# ------------------------------------------------------------------------
#%%

for frequency, timevar in {'monthly': 'DATE', 'interview': 'INTDATE', 'cuid': ''}.items():

    # levels of aggregation:
    agg_by = list(filter(None, ['CUID', timevar]))

    # aggregates rebates to monthly
    dfagg = aggregate_df(df.reset_index(), 
                               agg_by = agg_by, 
                               sum = sum_vars, 
                               max = max_vars,
                               min = min_vars)

    # merge with all household identifiers
    dfagg = dfagg.merge(hhdates.groupby(agg_by).max(), how='right', left_index=True, right_index=True)

    # if we are in the monthly or interview dataframe
    if timevar:
        # the time the last rebate was received (in monthly or interview) or the first date
        dfagg['LAST RBT ' + timevar] = dfagg['RBT ' + timevar]
        dfagg = dfagg.drop(min_vars, axis=1)

        # forward fill last dates
        dfagg['LAST RBT ' + timevar] = dfagg['LAST RBT ' + timevar].groupby('CUID').fillna(method='ffill')

        # fill missing values
        dfagg = fill_missing_by_type(dfagg)

        # Create Leads and Lags of Rebate Set to zero for March/April 2008, missing otherwise
        if frequency =='monthly':
            timeshift = 1           # expenditures are one month apart
        elif frequency =='interview':   
            timeshift = 3           # interviews are three months apart

        lagnum = int(6 / timeshift)
        leadnum = int(3 / timeshift)           

        # these are the variables we create lags of
        laglist = ['RBT INDICATOR','RBTAMT']
        dfshift = dfagg[laglist].reset_index()

        for i in range(-leadnum,lagnum+1):

            if i==0:
                continue
            elif i<0:
                prefix = 'LEAD'
            elif i>0:
                prefix = 'LAG'

            dflag = dfshift.copy()
            dflag[timevar] = dflag[timevar] + pd.DateOffset(months=i*timeshift)
            dflag = dflag.set_index(['CUID',timevar]).add_prefix(prefix + str(abs(i)))

            # merge with data
            dfagg = dfagg.join(dflag)

            # Replace nans with zeros if time < April 2008 or time > March 2009
            # this is the date of the lag or lead
            lagdate = dfagg.index.get_level_values(timevar) + pd.DateOffset(months= -i * timeshift)

            for var in laglist:
                if dfshift.dtypes[var] == 'bool':
                    dfagg.loc[(lagdate < '2008-04-01'), prefix + str(abs(i)) + var]   = False 
                    dfagg.loc[(lagdate > '2009-03-01'), prefix + str(abs(i)) + var]   = False
                else:
                    dfagg.loc[(lagdate < '2008-04-01'), prefix + str(abs(i)) + var]   = 0 
                    dfagg.loc[(lagdate > '2009-03-01'), prefix + str(abs(i)) + var]   = 0 

        # testing rebate construction
        if frequency =='interview': 
            assert dfagg.loc[pd.IndexSlice[186068,'2008-09-01'],'LAG1RBT INDICATOR'].all()
            assert dfagg.loc[pd.IndexSlice[186068,'2008-03-01'],'LEAD1RBT INDICATOR'].all()
            assert dfagg.loc[pd.IndexSlice[186068,'2007-12-01'],'LAG1RBT INDICATOR'].all()==False
        if frequency =='monthly': 
            assert dfagg.loc[pd.IndexSlice[186068,'2008-06-01'],'LAG1RBT INDICATOR'].all()
            assert dfagg.loc[pd.IndexSlice[186068,'2008-04-01'],'LEAD1RBT INDICATOR'].all()
            assert dfagg.loc[pd.IndexSlice[186068,'2008-02-01'],'LAG1RBT INDICATOR'].all()==False   
       
        
                     
    # if we are in the household dataframe
    else:
        # since we take the min, these are the first rebate date
        for var in min_vars:
            dfagg = dfagg.rename(columns={var: 'FIRST ' + var})
    
        # fill missing values
        dfagg = fill_missing_by_type(dfagg)

        # for household totals rename the variables: totals for sums, ever for maxes
        for var in dfagg.columns:
            if dfagg[var].dtype in ['int64', 'float64']:
                dfagg = dfagg.rename(columns={var: 'TOTAL ' + var}) 

            elif dfagg[var].dtype in ['bool', 'object']:     
                dfagg = dfagg.rename(columns={var: 'EVER ' + var}) 

        #Creates Categorical variable (0) Non-recipient, (1) Check only recipient, (2) EFT only recipient, (3) Received both, (9) RBT recipient but missing EFT/Check info
        
        dfagg['cat_eft'] = dfagg['EVER EFT'].multiply(1) + dfagg['EVER CHECK'].multiply(1) + 1
        dfagg.loc[((dfagg['EVER EFT']==False) & (dfagg['EVER CHECK']==True )), 'cat_eft'] = 1 
        dfagg.loc[dfagg['EVER RBT INDICATOR']==False, 'cat_eft'] = 0 
        dfagg.loc[((dfagg['EVER RBT INDICATOR']==True) & (dfagg['EVER EFT']==False) & (dfagg['EVER CHECK']==False )), 'cat_eft'] = 9

    # save file as parquet
    dfagg.to_parquet('../output/rebate' + frequency + '.parquet')
    
    if frequency =='monthly':
        dfsmall = dfagg[['RBTAMT','REBATES']]
        dfsmall.to_stata('../output/rebate' + frequency + '.dta')


#%%
