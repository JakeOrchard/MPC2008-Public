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

#import_from_nipa =  ['mv_parts', 'net_used', 'cars_nu','cars_n','cars_downpayment','pcars_n','trucks_n', 'usedautos','usedlighttrucks', 'home_dur', 'rec_dur',
#                     'other_dur', 'durable_nipa', 'foodbev_home',
#    Now Imports                  'clothing', 'gasoline', 'other_nondur', 'nondurable_nipa', 
#       all columns               'house_nipa', 'health_nipa', 'transportation', 'rec_service',
#                       'acc_service', 'fin_service', 'other_service', 'service_nipa','pce','tuition','mv_parts_ser', 'sndexp_nipa']


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
    df = pd.read_parquet('../output/psmjexpenditures' + expenditure_frequency + '.parquet')

    # load rebate information that varies at household level
    rebate = pd.read_parquet('../input/rebate' + expenditure_frequency + '.parquet')

    #Load NIPA Variables and convert index to uppercase
    nipa =  pd.read_parquet('../input/nipavariables' + expenditure_frequency + '.parquet')                    
    nipa.index = nipa.index.rename([x.upper() for x in list(nipa.index.names)])

   
    # merge family characteristics and rebate information
    for file_to_merge in interview, cohort, rebate_cohort, rebate,  nipa:
        df = df.merge(file_to_merge, how='left', left_index=True, right_index=True)



    # ------------------------------------------------------------------------
    # Save output
    # ------------------------------------------------------------------------

    df.to_parquet('../output/psmjvariables' + expenditure_frequency + '.parquet')
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
