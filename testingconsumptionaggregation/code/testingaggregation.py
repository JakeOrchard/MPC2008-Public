import pandas as pd
import yaml 

# ------------------------------------------------------------------------
# LOAD MTBI INPUT FILE
# ------------------------------------------------------------------------

# read aggregated MTBI file
dfmtbi = pd.read_parquet('../input/mtbimonthly.parquet')

# aggregate at interview level to compare with fmli files
dfaggall = dfmtbi.groupby(['NEWID']).sum()

# get expenditure names to compare
with open('../input/cexname_conversion.yml', 'r') as yamlfile:
        cexname_conversion = yaml.load(yamlfile, Loader=yaml.FullLoader)

# ------------------------------------------------------------------------
# LOAD FMLI INPUT FILE
# ------------------------------------------------------------------------

# loads fmli file from input folder
# parquet understands indexing so reindexing is not necessary
dffmli = pd.read_parquet('../input/fmliquarterly.parquet')

# ------------------------------------------------------------------------
# MERGE FRAMES AND CHECK
# ------------------------------------------------------------------------

dfmerge = dfaggall.merge(dffmli, left_on='NEWID', right_on='NEWID', suffixes=(None, '_F'))

# we write all the diagnostics to text as well
with open('../output/testingaggregation.txt', 'w') as f:

    for key in cexname_conversion:
        
        # differences in values
        dfmerge['Test'] =   (       
                                    dfmerge[cexname_conversion[key]]
                                -   dfmerge[cexname_conversion[key] + '_F']
                            )

        f.write(str(cexname_conversion[key]) + '\n')  
        f.write(str(dfmerge['Test'].describe()) + '\n')   

        print(dfmerge[abs(dfmerge['Test'])>1].shape[0])
        print(cexname_conversion[key])                 
        print(dfmerge['Test'].describe())
