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

filetype = ['nipa','psmj']

for file in filetype:
    
    if file == 'psmj':
        uccmap = '../input/ucc_category_map.yml'
        aggtype = '../input/psmjcategories_map.yml'
        totexp = 'TOTEXP2'

    elif file == 'nipa':
        uccmap = '../input/ucc_nipa_map.yml'
        aggtype = '../input/nipacategories_map.yml'
        totexp = 'pce'



    

    # ------------------------------------------------------------------------
    # Loads levels of PSMJ and NIPA variables
    # 1 - Total Expenditure
    # 2 - Major Subcategories such as NDEXP in PSMJ
    # 2.5 - Non-mutally exclusive major subcategories
    # 3 - Large subcategories such as New and Used cars in NIPA 
    # 4 - Lowest level of aggregation above pure UCC codes
    # ------------------------------------------------------------------------

    # Level 4 Labels
    with open(uccmap, 'r') as yamlfile:
        uccmapping = yaml.load(yamlfile, Loader=yaml.FullLoader)

    df= pd.DataFrame.from_dict(uccmapping.keys())
    df =df.rename(columns={0: "coefname"})    
    df['cattype'] = file
    df['agg_level'] = 4

    #Level 3 Labels PSMJ
    if file == 'psmj':
        with open('../input/level3to2_map.yml', 'r') as yamlfile:
            ucc32 = yaml.load(yamlfile, Loader=yaml.FullLoader)

        df32= pd.DataFrame.from_dict(ucc32.keys())
        df32 =df32.rename(columns={0: "coefname"})    
        df32['cattype'] = file
        df32['agg_level'] = 3

       


    #Level 2 Labels PSMJ/ Level 3 Labels in NIPA
    with open(aggtype, 'r') as yamlfile:
        aggmapping = yaml.load(yamlfile, Loader=yaml.FullLoader)

    dfagg = pd.DataFrame.from_dict(aggmapping.keys())
    dfagg =dfagg.rename(columns={0: "coefname"})    
    dfagg['cattype'] = file
    dfagg['agg_level'] = 2
    if file == 'nipa':
        dfagg['agg_level'] = 3
        #Separate dataframe dictionary that maps level 4 to level 3, and from level 3 to 2
        df4to3 = pd.DataFrame.from_dict(aggmapping, orient='index').melt(ignore_index = False , value_name='Level4')
        df4to3 = df4to3[df4to3['Level4'].notnull()]
        df3to2 = df4to3.loc[['NONDURABLE_NIPA','DURABLE_NIPA', 'SERVICE_NIPA']].reset_index()
        df3to2 = df3to2.rename(columns={"index": "Level2", "Level4": "Level3"})
        #PULLS OUT PSMJ NONDURALBES In PCE SPACE
        dfnd = df4to3.loc[['NDEXP_NIPA','SNDEXP_NIPA']].reset_index()
        dfnd = dfnd.rename(columns={"index": "ND_PSMJ"})
        #Drops levels 1 and 2 from 4to3 dictionary, and merges in 3 to 2 dictionary
        df4to3 = df4to3[~df4to3.index.isin(['NONDURABLE_NIPA','DURABLE_NIPA', 'SERVICE_NIPA','PCE'])].reset_index()
        df4to3 = df4to3.rename(columns={"index": "Level3"})
        df4to2 = df3to2.merge(df4to3,how='inner',on='Level3')
        df4to2 = df4to2.merge(dfnd,how='outer',on='Level4')
        df4to2 = df4to2.loc[:, ~df4to2.columns.str.startswith('variable')]
        df4to2['coefname'] = df4to2['Level4']
        df4to2.to_stata('../output/agglabels_levelnames_' + file + '.dta')



#%%
   #Append together
    df = pd.concat([df,dfagg])

   #PSMJ renames some variables, and corrects for non-mutually exclusive categories
    if file == 'psmj':
        df = pd.concat([df,df32])
        df = df.set_index('coefname')
       
        with open('../input/cexname_conversion.yml', 'r') as yamlfile:
                cexname_conversion = yaml.load(yamlfile, Loader=yaml.FullLoader)
        
        df = df.rename(index=cexname_conversion)
        df.loc[['FOODBEVS','SNDEXP','USED_CAR_PCE','TOTEXP_NMV'],'agg_level'] = 2.5  #Level 2's are not mutually exclusive in PSMJ space, relabelling some as level 2.5
        df.loc[['VEHINS','FDHOME', 'FDAWAY', 'TOBACC', 'PUBTRA', 'GASMO', 'READ', 'CARTKN', 'CARTKU', 'OTHVEH', 'HOUSEQ', 'VEHFIN', 'EDUCA', 'MAINRP', 'VEHFIN'], 'agg_level'] = 3.5
        df.loc[['PHONED',  'CASHCONT', 'NEWPCARS', 'usedlighttrucks', 'BOOKS', 'DRUGS', 'MAGAZINES', 'RETPEN', 'usedautos', 'VEHPURCH', 'NEWTRUCKS',  'LIFINS'], 'agg_level'] = 4.5

    #Relabels Major Categories in NIPA Space
    if file == 'nipa':
        droprows = ['FOODHOME', 'FOODAWAY','ALCHOME','ALCAWAY']
        df = df[df.coefname.isin(droprows) == False]
        df = df.set_index('coefname')
        df.index = df.index.str.lower()
        df.loc[['durable_nipa','nondurable_nipa','service_nipa'], 'agg_level'] = 2  #Major Nipa categories
        df.loc[['foodbevs','ndexp_nipa','sndexp_nipa','mv_parts_ser','tuition','alcbever','trucks_n','pcars_n','cars_downpayment','net_used','cars_nu','cars_n'], 'agg_level'] = 2.5  #Major Nipa categories
        df.loc[['newcardown','usedcardown','mv_parts_ser','tuition','alcbever','trucks_n','pcars_n','cars_downpayment','net_used','cars_nu','cars_n'], 'agg_level'] = 2.5  #Major Nipa categories
        df.loc[[], 'agg_level'] = 4.5


    #Level 1 Label
    df.loc[[totexp], 'agg_level'] = 1 #
    
    #Saves as Stata file
    df = df.reset_index()
    df['coefname'] = df['coefname'].str.lower()
    df.to_stata('../output/agglabels' + file + '.dta')
    #%%
# %%
