#%%
import zipfile
import os, fnmatch
import pandas as pd
import numpy as np
import shutil

#%%
# extract zip file to this folder
zipextractdir = '../output/cexdownload/'

# directory of downloaded files
zipfilepath = '../input/cexdownload.zip'

# unzip files
with zipfile.ZipFile(zipfilepath, 'r') as zip_ref:
        zip_ref.extractall(zipextractdir)

#%%
# function that finds all the files that match a certain pattern
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

#%% 
# for each file type we find the files starting with the pattern in square 
# brackets 
for filetype in ["fmli", "mtbi", "rbt"]:

    # find files of the type that end with 2,3,4, or x
    # this excludes the first quarter files without x
    if filetype=="fmli" or filetype=="mtbi" :
        files_to_append = find(filetype + '*[234x].csv', zipextractdir) 
    else:
        files_to_append = find(filetype + '*.csv', zipextractdir) 
    
    print(files_to_append)

    # this code appends the files in a pandas dataframe
    df = pd.concat(map(pd.read_csv, files_to_append))

    # need to fix some formatting issues in the FMLI files
    if filetype=='fmli':
        # this row is incomplete and misaligned in the raw data:
        #df.drop(df.loc[df['NEWID']==2361325].index, inplace=True)

        # replace "." entries with empty strings to allow conversion to numeric
        df = df.replace(to_replace='.', value='')
        df['PSU'] = df['PSU'].replace(to_replace=['S[0-9][0-9][A-Z]'], value=np.nan, regex=True)

    # applies numeric conversion when feasible
    df = df.apply(pd.to_numeric, errors='ignore')

    # save as parquet file
    df.to_parquet('../output/' + filetype + '.parquet')

    # print(df)

# %%
# delete directory with extracted files    
shutil.rmtree(zipextractdir)

# %% Save as stata file

df.to_stata('../output/df' + '.dta')

# %%
