#%%
import pandas as pd

# ------------------------------------------------------------------------
# Loading and aggregating FMLI files
# ------------------------------------------------------------------------

fmli = pd.read_parquet('../input/fmliquarterly.parquet')
fmli['YEAR'] = pd.DatetimeIndex(fmli['INTDATE']).year



# ------------------------------------------------------------------------
# Create cohort variables (fixed at CUID level)
# ------------------------------------------------------------------------

# have to do this differently
# create a variable that captures first date in sample
firstdate = fmli.groupby('CUID')['INTDATE'].min().rename('FIRSTDATE')

# create a variable that captures whether a given CUID was present for
# interview 2,3,4,5
interviews = fmli[['CUID','INTNUM']].set_index(['CUID','INTNUM'])

interviews['EVER'] = 1


interviews = interviews.unstack().fillna(0)
print(interviews)
interviews.columns = [ 'EVER INT 1', 'EVER INT 2','EVER INT 3','EVER INT 4','EVER INT 5']


# merge characteristics:
cohort = pd.DataFrame(index = firstdate.index)
for df in firstdate, interviews:
    cohort = cohort.join(df)

#%%
# ------------------------------------------------------------------------
# Save Output
# ------------------------------------------------------------------------

cohort.to_parquet('../output/cohortvariables.parquet')

