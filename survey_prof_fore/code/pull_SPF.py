
# %%

#%% Packages
import pandas as pd
import numpy as np
from time import strptime



#%%  Download all sheets

url = 'https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/files/individual_rconsum.xlsx'


df = pd.read_excel(url) 




#%% Labels industries

df['INDUSTRY_LABEL'] = 'Financial Service Provider'
df.loc[df['INDUSTRY']== 2,'INDUSTRY_LABEL'] = 'Non-Financial Service Provider'
df.loc[df['INDUSTRY']== 3,'INDUSTRY_LABEL'] = "Unknown"
df.loc[df['INDUSTRY']== np.nan,'INDUSTRY_LABEL'] = "Unknown"

# %%

df.to_stata('../output/spf_rcons.dta')


# %%
