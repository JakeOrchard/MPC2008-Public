#%%
import pandas as pd
import numpy as np

# ------------------------------------------------------------------------
# Loading and aggregating MTBI files
# ------------------------------------------------------------------------

# expenditure variables to be used in regression
# for the first set we will also add the log of expenditure
keep_variables_dur = ['MV_PARTS','NET_USED','CARS_NU','CARS_N', 'CARS_DOWNPAYMENT','PCARS_N','TRUCKS_N','usedautos','usedlighttrucks','HOME_DUR','REC_DUR','OTHER_DUR','DURABLE_NIPA']
keep_variables_nondur = ['FOODBEV_HOME','CLOTHING','GASOLINE','OTHER_NONDUR','NONDURABLE_NIPA']                 
keep_variables_service = ['HOUSE_NIPA','HEALTH_NIPA','TRANSPORTATION','REC_SERVICE','ACC_SERVICE', 'FIN_SERVICE', 'OTHER_SERVICE','SERVICE_NIPA','PCE']
keep_variables_extra = ['TUITION','MV_PARTS_SER', 'SNDEXP_NIPA']
keep_variables_all = keep_variables_dur+keep_variables_nondur+keep_variables_service + keep_variables_extra

# read in monthly expenditure files all columns
dfmtbi = pd.read_parquet('../input/mtbimonthly_NIPACAT.parquet')
dfmtbi = dfmtbi.drop(columns=['CUID', 'index', 'INTCOUNT','']) 
#Remove food columns used for sample selection in create consumption variables
removelist = ['FOODHOME','FOODAWAY','ALCHOME','ALCAWAY']
dfmtbi = dfmtbi.drop(removelist, axis=1)
dfmtbi = dfmtbi.add_prefix('nipa_')

# add interview date to dataframe and set index to household - interview date
dffmli = pd.read_parquet('../input/fmliquarterly.parquet', columns = ['NEWID','CUID','INTDATE'])




#%%
#Merge Expenditure Data with family characteristics
dfmtbi = dfmtbi.merge(dffmli, how='left', right_index=True, left_index=True)


#%%
# ------------------------------------------------------------------------
# Loop over aggregation and save data
# ------------------------------------------------------------------------


for expenditure_frequency in ['monthly', 'interview']: 

    # existing data already at monthly frequency
    if expenditure_frequency=='monthly':
        timevar = 'DATE'          
        
        # set index to household ID and expenditure month
        df = dfmtbi.reset_index().set_index(['CUID', 'NEWID', timevar, 'INTDATE']) 


    # for interview level we aggregate and merge interview date
    elif expenditure_frequency=='interview':   
        timevar = 'INTDATE'

        # aggregate to interview level 
        df = dfmtbi.groupby(['CUID', 'NEWID', timevar]).sum()


   
    df.to_parquet('../output/nipaexpenditures' + expenditure_frequency + '.parquet')

#%%
# ------------------------------------------------------------------------
# Testing aggregation
# ------------------------------------------------------------------------
# if __name__=='__main__':
#     df = pd.read_parquet('../output/nipaexpendituresmonthly.parquet', columns = keep_variables_all)
#     df = df.groupby(['CUID', 'NEWID']).sum()    

#     dfint = pd.read_parquet('../output/nipaexpendituresinterview.parquet', columns = keep_variables_all)
#     dfint = dfint.add_suffix('_int')
#     # df = dfint.groupby(['CUID', 'NEWID']).sum()

#     df = df.merge(dfint, how='left', left_index=True, right_index=True)

#     for var in keep_variables_all:
#         print('Aggregation Residuals for ' + var + ':')
#         print((df[var]-df[var + '_int']).abs().sum())
#         print((df[var].isna()^df[var + '_int'].isna()).abs().sum())


# %%
