
# %%

#%% Packages
import pandas as pd
import numpy as np
from time import strptime
import wget
import os
import zipfile


#%%  Download and unzip forecast

filename = "greenbook_rpce"

url = 'https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/greenbook-data/gbweb/gbweb_grpce_column_format.zip?la=en&hash=15BCE1E604E3B5EFA955E279D751EB45'


# construct directory
zipfilepath = '../output/'  + filename 

# extract zip file to this folder
zipextractdir = '../output/greenbookforecast/'

wget.download(url, zipfilepath)


with zipfile.ZipFile(zipfilepath, 'r') as zip_ref:
    zip_ref.extractall(zipextractdir)

    # delete downloaded zip file
    os.remove(zipfilepath)


#%%
excelfilename = zipextractdir + 'gRPCE_1985_Last.xlsx'
df = pd.read_excel(excelfilename) 







# %%

df.to_stata('../output/greenbook_rpce.dta')


# %%
