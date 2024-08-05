#%%
import pandas as pd

# ------------------------------------------------------------------------
# Variables to Import
# ------------------------------------------------------------------------

# household characteristics that vary by interview
import_from_interview = ['AGE','d_AGE',
                        'NUM_ADULTS','PERSLT18','FAM_SIZE','d_NUM_ADULTS','d_PERSLT18',
                        'CUTENURE',
                        'FINLWT21','FINCBTAX','LINCOME','INCD INTERVIEW', 'INCD POOL']

interview = pd.read_parquet('../input/interviewvariables.parquet', 
                            columns = import_from_interview)                    

# variables fixed by household
cohort = pd.read_parquet('../input/cohortvariables.parquet') 

# rebate characteristics fixed by household
rebate_cohort = pd.read_parquet('../input/rebatecuid.parquet')
#%%
# ------------------------------------------------------------------------
# Merge each level of aggregation
# ------------------------------------------------------------------------

# could run two versions, monthly or interview level
for expenditure_frequency in ['monthly', 'interview']:  #

    # read expenditure data (already logged and in difference, varying by frequency)
    df = pd.read_parquet('../output/nipaexpenditures' + expenditure_frequency + '.parquet')

    # load rebate information that varies at household level
    rebate = pd.read_parquet('../input/rebate' + expenditure_frequency + '.parquet')

    # merge family characteristics and rebate information
    for file_to_merge in interview, cohort, rebate_cohort, rebate:
        df = df.merge(file_to_merge, how='left', left_index=True, right_index=True)



    # ------------------------------------------------------------------------
    # Save Output
    # ------------------------------------------------------------------------



    # save combined file as csv to output folder and then zip
    df.to_parquet('../output/nipavariables' + expenditure_frequency + '.parquet')

    # converts all True/False booleans to numeric
    for col in df.select_dtypes(include='bool').columns:
       df[col] = df[col].astype(int)

    for col in df.select_dtypes(include='object').columns:
       df[col] = df[col].multiply(1)
       df[col] = df[col].astype('float64')

   # make all columns and index names lower case
    df.columns = [col.lower().replace(' ', '') for col in df.columns ]
    df.index = df.index.rename([x.lower() for x in list(df.index.names)])

   # save as stata and csv
    dates_format = {'intdate': 'tm','firstrbtdate':'tm','firstrbtintdate':'tm'}
    if expenditure_frequency=='monthly':
       dates_format['date'] = 'tm'

    df.to_stata('../output/nipavariables' + expenditure_frequency + '.dta',
                  convert_dates = dates_format)

    df.to_parquet('../output/nipavariables' + expenditure_frequency + '.parquet')

    
#%%

# testing whether file is the same as our old version
# fixed some errors that were present in the old version, so no longer informative
# if __name__=='__main__':
#     df = pd.read_parquet('../output/psmjvariablesinterview.parquet')
#     dfold = pd.read_parquet('../output/psmjvariables.parquet')


#     for col in df.columns:
#         print(col)
#         if col in dfold.columns:
#             print((df[col]-dfold[col]).abs().sum())
#             print((df[col].isna()^dfold[col].isna()).abs().sum())
#             print(df[col].isna().abs().sum())
#             print(dfold[col].isna().abs().sum())

# %%
