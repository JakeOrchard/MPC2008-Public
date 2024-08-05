#%%
import pandas as pd
import yaml 

import os, sys
rootpath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(rootpath)
from utilities.mpc_utilities import aggregate_with_dict

# ------------------------------------------------------------------------
# Imports CE PUMD dictionary to name missing UCCs
# ------------------------------------------------------------------------

uccdict = pd.read_excel("../input/ce-pumd-interview-diary-dictionary.xlsx",
                        sheet_name = 'Codes ',
                        header = 0,
                        index_col = [0,1,2])

uccdict = uccdict.loc[pd.IndexSlice['INTERVIEW','MTBI','UCC'], :]
uccdict = uccdict.reset_index()
uccdict['Code value'] = pd.to_numeric(uccdict['Code value']) 
uccdict = uccdict[(uccdict['Last year']<1996) == False]
uccdict = uccdict[(uccdict['Code description']=='VEHICLE REGISTRATION STATE/LOCAL') == False]
uccdict = uccdict.set_index('Code value')
cecodes = uccdict.index.get_level_values(0)

# ------------------------------------------------------------------------
# Loop create variables for PSMJ Replication and then PCE-NIPA comparison
# ------------------------------------------------------------------------

filetype = ['nipa','psmj']

for file in filetype:
    
    if file == 'psmj':
        read_file = '../input/mtbi.parquet'
        uccmap = '../input/ucc_category_map.yml'
        aggtype = '../input/psmjcategories_map.yml'
        savefile = '../output/mtbimonthly.parquet'
        totexp = 'TOTEXP2'

    elif file == 'nipa':
        read_file = '../output/mtbi_wnewUCC.parquet'
        uccmap = '../input/ucc_nipa_map.yml'
        aggtype = '../input/nipacategories_map.yml'
        savefile = '../output/mtbimonthly_NIPACAT.parquet'
        totexp = 'PCE'



    
    # ------------------------------------------------------------------------
    # LOAD INPUT FILE
    # ------------------------------------------------------------------------

    # reads mtbi parquet file
    df = pd.read_parquet(read_file)

    
    # create expenditure date
    df['DATE'] = pd.to_datetime(dict(year = df['REF_YR'], 
                                    month = df['REF_MO'], 
                                    day = 1))
    
    # drop columns we do not use
    df.drop(['COST_','GIFT','PUBFLAG','REF_MO','REF_YR','SEQNO','ALCNO','UCCSEQ','EXPNAME','RTYPE'], inplace=True, axis=1)

    # ------------------------------------------------------------------------
    # COLLAPSE EXPENDITURE BY ID, UCC, AND DATE
    # ------------------------------------------------------------------------

    # aggregate UCCs
    df = df.groupby(['NEWID','DATE','UCC']).sum() 


    # ------------------------------------------------------------------------
    # DATA IS AT UCC CODE LEVEL. THE NEXT STEP CREATES THE MAPPING FROM UCC TO
    # CATEGORIES AND THEN SUMS TO MONTHLY
    # ------------------------------------------------------------------------

    # import YAML file with ucc mapping to category
    with open(uccmap, 'r') as yamlfile:
        uccmapping = yaml.load(yamlfile, Loader=yaml.FullLoader)

    # loop over categories
    df['CAT'] = ''
    uccsindata = df.index.get_level_values('UCC').unique()
    uccsinmapping = []
    for catcode, ucclist in uccmapping.items():
        for ucc in ucclist:
            uccsinmapping.append(ucc)
    with open('../output/' + file + '_ucc_not_mapped.txt', 'w') as f:
        for ucc in uccsindata:
            if ucc not in uccsinmapping and ucc<790000 and ucc>10000:
                if ucc in cecodes:
                    f.write(str(ucc) + ' ' + uccdict._get_value(ucc, 'Code description')+ ', ' + str(uccdict._get_value(ucc, 'First year')) + ' ' + str(uccdict._get_value(ucc, 'Last year')) + '\n') 
                else:
                    f.write(str(ucc) + '\n') 

    for catcode, ucclist in uccmapping.items():
        # for catcode, ucclist in category.items():
        # need to check if ucc has already been assigned 
        print(catcode)
        for ucc in ucclist:
            if ucc not in uccsindata:
                continue
            if (df.loc[pd.IndexSlice[:,:,ucc],['CAT']] == '').all()[0]:
                # this creates NIPA entries whenever the ucc code falls in the right 
                df.loc[pd.IndexSlice[:,:,ucc],['CAT']] = catcode
            else:
                print(ucc)
                print('making copy')
                # if already assigned, create duplicates of uccs and assign new code
                dfslicecopy = df.loc[pd.IndexSlice[:,:,ucc],:].copy()
                dfslicecopy['CAT'] = catcode
                df = df.append(dfslicecopy)


    # for category in uccmapping:

    #     # this creates cost entries whenever the ucc code falls in the right 
    #     # category
    #     df[category] = df['COST'] * df['UCC'].isin(uccmapping[category])

    # drop columns we do not use
    df = df.reset_index().drop(['UCC'], axis=1)

    # aggregate at monthly frequency
    df = df.fillna(0).groupby(['NEWID','DATE','CAT']).sum()['COST']
    df = df.unstack().fillna(0)
    for catcode in uccmapping.keys():
        if catcode not in df.columns:
            df[catcode] = 0
    
    # ------------------------------------------------------------------------
    # CREATE CONSUMPTION AGGREGATES
    # ------------------------------------------------------------------------



    if file == 'psmj':
        # create level two consumption expenditure categories
        df = aggregate_with_dict(df,'../input/level3to2_map.yml') 

        # rename expenditure categories to short
        with open('../input/cexname_conversion.yml', 'r') as yamlfile:
                cexname_conversion = yaml.load(yamlfile, Loader=yaml.FullLoader)

        df.rename(columns=cexname_conversion, inplace=True)

 
    # aggregate to expenditure variables
    df = aggregate_with_dict(df,aggtype) 

    
    # ------------------------------------------------------------------------
    # CLEAN AGGREGATED DATA
    # ------------------------------------------------------------------------

    if file == 'psmj':
        foodvar = 'FOODBEVS'
    elif file == 'nipa':
        foodvar = 'FOODBEVS'

    # drop rows with a fourth or fifth interview and zero food expenditure and all rows with no expenditure
    df = df.merge(df[foodvar].groupby('NEWID').count().rename('INTCOUNT'), how='left', left_index=True, right_index=True)
    df = df[((df[foodvar].abs()>0) | (df['INTCOUNT']<=3)) & (df[totexp].abs()>0)]
    
    #Checking to see if NIPA drops too many rows (rows where households have no at home food expenditure)

    # check that the index is unique with month and CUID:
   
    dftest = df[['ALCAWAY','ALCHOME']]

    dftest = dftest.reset_index()

    # CUID value is NEWID dropping the last digit (which is interview number)
    dftest['CUID'] = dftest['NEWID'] // 10

    dftest = dftest.set_index(['CUID','DATE'])

    
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


    
    # ------------------------------------------------------------------------
    # SAVE OUTPUT
    # ------------------------------------------------------------------------
    # checks that panel is unique
   
    assert dftest.index.value_counts().max()==1, 'Expenditure file does not have unique panel identifiers. Do not continue.'



# ------------------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------------------

    # save combined file to parquet
    df.to_parquet(savefile, index=True)


# %%
