#%% Packages
import pandas as pd
import numpy as np
from time import strptime



#%%  Download all sheets

url = 'https://www.census.gov/retail/mrts/www/mrtssales92-present.xls'

headersize = 4
footersize = 9
df = pd.read_excel(url, sheet_name = None,skiprows= range(0, headersize)) # read all sheets


# %% Only keep adjusted series and transpose


#Most recent year is in different format. 
max_year = int(list(df)[0])

#Only keep SA adjusted series 
df_adjusted = df[str(max_year)][67:105] #SA series starts here for each. Can generalize later
df_adjusted= df_adjusted.drop(columns=['Unnamed: 0', 'CY CUM', 'PY CUM'])
df_adjusted.rename({'Unnamed: 1': 'Category'}, axis=1, inplace=True)
df_adjusted = df_adjusted.set_index('Category')

#Reshape data from wide to long
df_adjusted = df_adjusted.stack().reset_index()
df_adjusted.rename({'level_1': 'month',0:'sales'}, axis=1, inplace=True)

#Will place cleaned data in this frame
df_clean = df_adjusted

#%% Now loops through the rest of the years, which have common format

# Current year (2022) has slightly different format
for i in range(1992,max_year):
    
    #Only keep SA adjusted series 
    df_adjusted = df[str(i)][67:] #SA series starts here for each. Can generalize later
    df_adjusted= df_adjusted.drop(columns=['Unnamed: 0', 'TOTAL'])
    df_adjusted.rename({'Unnamed: 1': 'Category'}, axis=1, inplace=True)
    df_adjusted = df_adjusted.set_index('Category')
    
    #Reshape data from wide to long
    df_adjusted = df_adjusted.stack().reset_index()
    df_adjusted.rename({'level_1': 'month',0:'sales'}, axis=1, inplace=True)
    df_clean = df_clean.append(df_adjusted)


#%%
#Convert months to date-time
df_clean[['m','year']]=df_clean['month'].str.split(r"\s", n=-1,expand=True)[[0,1]]
maxstring = str(max_year) 
maxstringp = maxstring + "(p)"
df_clean['year'][df_clean['year'] == maxstringp] = maxstring 
df_clean['m'] = df_clean['m'].str.slice(stop=3)
df_clean['month'] = df_clean['m'] + df_clean['year']
df_clean['date'] = pd.to_datetime(df_clean['month'], infer_datetime_format=True)
df_clean = df_clean.drop(columns=['year', 'm','month'])
df_clean = df_clean.replace(to_replace="(S)",value = np.nan)
df_clean.sort_values(by=['Category', 'date'])
df_clean = df_clean.set_index(['date','Category'])
df_clean = df_clean.unstack()

#%% Convert Months to python time

#%%

#%% Export data to csv and stata

df_clean.to_csv('../output/retail_sales_clean.csv')

#%%
#dates_format = {'date': 'tm'}
df_clean.to_stata('../output/retail_sales_clean.dta')

    





# %%
