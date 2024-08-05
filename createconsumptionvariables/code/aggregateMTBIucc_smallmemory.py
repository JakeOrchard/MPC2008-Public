#%%
import pandas as pd
import yaml 

import os, sys
rootpath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(rootpath)
from utilities.mpc_utilities import aggregate_with_dict

# ------------------------------------------------------------------------
# Loop create variables for PSMJ Replication and then PCE-NIPA comparison
# ------------------------------------------------------------------------

filetype = ['psmj','nipa']

for type in filetype:
    for year in range(1996,2020):
        if type == 'psmj':
            read_file = '../input/mtbi.parquet'
            uccmap = '../input/ucc_category_map.yml'
            aggtype = '../input/psmjcategories_map.yml'
            savefile = '../output/mtbimonthly' + str(year) + '.parquet'
            totexp = 'TOTEXP2'

        elif type == 'nipa':
            read_file = '../output/mtbi_wnewUCC.parquet'
            uccmap = '../input/ucc_nipa_map.yml'
            aggtype = '../input/nipacategories_map.yml'
            savefile = '../output/mtbimonthly_NIPACAT' + str(year) + '.parquet'
            totexp = 'PCE'



        #%%
        # ------------------------------------------------------------------------
        # LOAD INPUT FILE
        # ------------------------------------------------------------------------

        # reads mtbi parquet file
        df = pd.read_parquet(read_file)
        df = df[df['REF_YR']==year]


        #%% ------------------------------------------------------------------------
        # DATA IS AT UCC CODE LEVEL. THE NEXT STEP CREATES THE MAPPING FROM UCC TO
        # CATEGORIES AND THEN SUMS TO MONTHLY
        # ------------------------------------------------------------------------

        # import YAML file with ucc mapping to category
        with open(uccmap, 'r') as yamlfile:
            uccmapping = yaml.load(yamlfile, Loader=yaml.FullLoader)

        # loop over categories
        for category in uccmapping:

            # this creates cost entries whenever the ucc code falls in the right 
            # category
            df[category] = df['COST'] * df['UCC'].isin(uccmapping[category])

        # create expenditure date
        df['DATE'] = pd.to_datetime(dict(year = df['REF_YR'], 
                                        month = df['REF_MO'], 
                                        day = 1))


        # drop columns we do not use
        df.drop(['UCC','COST','COST_','GIFT','PUBFLAG','REF_MO','REF_YR'], inplace=True, axis=1)

        # aggregate at monthly frequency
        df = df.groupby(['NEWID','DATE']).sum()

        # ------------------------------------------------------------------------
        # CREATE CONSUMPTION AGGREGATES
        # ------------------------------------------------------------------------



        if type == 'psmj':
            # create level two consumption expenditure categories
            df = aggregate_with_dict(df,'../input/level3to2_map.yml') 

            # rename expenditure categories to short
            with open('../input/cexname_conversion.yml', 'r') as yamlfile:
                    cexname_conversion = yaml.load(yamlfile, Loader=yaml.FullLoader)

            df.rename(columns=cexname_conversion, inplace=True)

        #%%
        # aggregate to expenditure variables
        df = aggregate_with_dict(df,aggtype) 

        #%%
        # ------------------------------------------------------------------------
        # CLEAN AGGREGATED DATA
        # ------------------------------------------------------------------------

        if type == 'psmj':
            foodvar = 'FOODBEVS'
        elif type == 'nipa':
            foodvar = 'FOODBEV_HOME'

        # drop rows with a fourth or fifth interview and zero food expenditure and all rows with no expenditure
        df = df.merge(df[foodvar].groupby('NEWID').count().rename('INTCOUNT'), how='left', left_index=True, right_index=True)
        df = df[((df[foodvar].abs()>0) | (df['INTCOUNT']<=3)) & (df[totexp].abs()>0)]


        # check that the index is unique with month and CUID:
        if type == 'psmj':
            dftest = df[['ALCAWAY','ALCHOME']]
        elif type == 'nipa':
            dftest = df[['alcaway']]

        dftest = dftest.reset_index()

        # CUID value is NEWID dropping the last digit (which is interview number)
        dftest['CUID'] = dftest['NEWID'] // 10

        dftest = dftest.set_index(['CUID','DATE'])

        #%%
        # count the number of times NEWID shows up for each ID-Date combination
        # There appear to be 25 households that have repeated
        #CUID-DATE observations. None of these are in the rebate sample
        dftest['count'] = dftest.index.value_counts()

        if dftest['count'].max()>1:
            print('Expenditure file does not have unique panel identifiers.')

        #For now we drop these observations
        df = df.reset_index()
        df['CUID'] = df['NEWID']//10
        df.set_index(['CUID','DATE'])
        df['duplicates'] = df.index.value_counts()
        df = df[df['duplicates']==1]
        df = df.reset_index().set_index(['NEWID','DATE'])
        df= df.drop(columns={'duplicates'})


        #%%
        # ------------------------------------------------------------------------
        # SAVE OUTPUT
        # ------------------------------------------------------------------------
        # checks that panel is unique
        assert dftest.index.value_counts().max()==1, 'Expenditure file does not have unique panel identifiers. Do not continue.'



    #%% ------------------------------------------------------------------------
    # SAVE OUTPUT
    # ------------------------------------------------------------------------

        # save combined file to parquet
        df.to_parquet(savefile, index=True)
#%%
