# %%
import pandas as pd

# %%----------------------------------------------------------------------
# LOAD INPUT FILE
# ------------------------------------------------------------------------

# function that reads appended fmli files
df = pd.read_parquet('../input/fmli.parquet')

#df = df[(df['QINTRVYR']>=2007) &(df['QINTRVYR']<=2009)]

print('The Number of observations in the fmli file is (sumFMLIquarterly consumption.py):')
print(str(df.shape[0])) # There are 82909 in 2007-2009 sample and 82897 in the trimmed 1996-2019 sample

#All CUID values are NANs for pre-2003 data. Replace with NEWID/10

df['CUID'] = df['NEWID']//10



#%% ------------------------------------------------------------------------
# SET INDEX
# ------------------------------------------------------------------------ 
                   
# set index of dataframe
df.set_index('NEWID', inplace=True) 

# %%----------------------------------------------------------------------
# DATA IS SEPARATED BY PQ (PREVIOUS QUARTER) AND CQ (CURRENT QUARTER) BASED
# ON THE INTERVIEW MONTH. THE NEXT STEP SUMS THESE TWO VARIABLES TO
# CALCULATE TOTAL CONSUMPTION OVER PREVIOUS THREE MONTHS FROM INTERVIEW
# DATE
# ------------------------------------------------------------------------

# column names for the variables that end with CQ. these are the expenditure
# values for the previous quarter. 
df_columns_filter = df.filter(regex='CQ$',axis=1).columns

# next we sum previous quarter and current quarter expenditure variables
for column_name in df_columns_filter:

    # remove the CQ subscript from the column name
    expen_column = column_name.replace('CQ','',1)

    # SUM CQ and PQ expenditures
    df[expen_column] = df[expen_column + 'PQ'] + df[expen_column + 'CQ']

    # Drop CQ and PQ columns
    df.drop([expen_column + 'PQ',expen_column + 'CQ'], inplace=True, axis=1)

# %%----------------------------------------------------------------------
# CONSTRUCT INTERVIEW DATE AND INTERVIEW NUMBER
# ------------------------------------------------------------------------

# Interview date in which rebate is reported
df['INTDATE'] = pd.to_datetime(dict(year = df['QINTRVYR'], 
                                    month = df['QINTRVMO'], 
                                    day = 1))


df.drop(['QINTRVMO','QINTRVYR'], inplace=True, axis=1)  

# Interview number is last digit of NEWID, which is in the index
df['INTNUM'] = df.index % 10

# ------------------------------------------------------------------------
# INIDCATOR WHETHER REBATE MODULE WAS PART OF THE INTERVIEW
# ------------------------------------------------------------------------

df['RBTINTVIEW'] =  (
                      (df['INTDATE']>='2008-06-01') 
                    & (df['INTDATE']<='2009-03-01')
                    )



# %%----------------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------------------

# save to parquet file
df.to_parquet('../output/fmliquarterly.parquet', index=True)

