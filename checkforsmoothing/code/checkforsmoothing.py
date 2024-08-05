import pandas as pd
import yaml 

# ------------------------------------------------------------------------
# LOAD MTBI INPUT FILE
# ------------------------------------------------------------------------

# read aggregated MTBI file
dfmtbi = pd.read_parquet('../input/mtbimonthly.parquet')


# ------------------------------------------------------------------------
# WANT TO CHECK WHETHER EXPENDITURES WITHIN INTERVIEWS ARE MORE STRONGLY 
# CORRELATED THAN EXPENDITURES ACROSS INTERVIEWS
# ------------------------------------------------------------------------

# average expenditure in interview
dfavg = dfmtbi.groupby(['NEWID']).mean()

# checks whether any expenditure category is nonzero in ITW or ITW month:
dfitwnotzero = dfavg.ne(0).sum()
dfitwmonthnotzero = dfmtbi.ne(0).sum()

# checks whether expenditures are the same in any month
dfsame = ((dfmtbi-dfavg)/dfavg*100).abs().groupby(['NEWID']).mean().eq(0).sum()

# we write all the diagnostics to text as well
with open('../output/checkforsmoothing.txt', 'w') as f:

    for col in dfmtbi.columns:
        print(col)
        print('Nonzero Exp in ITW = ' + str(dfitwnotzero[col]))
        print('All Months the Same = ' + str(dfsame[col]))
        print('Nonzero Exp in ITW Month = ' + str(dfitwmonthnotzero[col]))
    
        f.write(str(col) + '\n')  
        f.write('Nonzero Expenditure in ITW = ' + str(dfitwnotzero[col]) + '\n') 
        f.write('All Months the Same in ITW = ' + str(dfsame[col]) + '\n')  
        f.write('Nonzero Expen in ITW Month = ' + str(dfitwmonthnotzero[col]) + '\n\n')