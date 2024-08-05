# %%

import pandas as pd
import requests
import json
import pandas_datareader.data as web    




# ------------------------------------------------------------------------
# Downloads Data from BEA
# ------------------------------------------------------------------------
UserID = 'INSERT BEA CODE HERE'
TableType = 'NIUnderlyingDetail'
TableNames = ['U20405','U20404','U20406']
TableTitles = ['pceexpenditure','pceexpenditureprice','pceexpenditurereal']
Frequency = 'M'
Yearstarts = [1959, 1959, 2007]
Yearend = 2022


for TableName, TableTitle, Yearstart in zip(TableNames, TableTitles, Yearstarts):
    # Loop over years such that JSON file is not too large for BEA API
    for Year in range(Yearstart,Yearend + 1):
        # print('Downloading Year ' + str(Year))

        url =   (
                'https://apps.bea.gov/api/data/?'
                + '&UserID=' + UserID 
                + '&method=GetData&'
                + 'DataSetName=' + TableType 
                + '&TableName=' + TableName 
                + '&Frequency=' + Frequency 
                + '&Year=' + str(Year) 
                + '&ResultFormat=json'
                )

        # beatable = requests.get(url, stream=False, proxies=proxies, headers=headers)
        beatable = requests.get(url)
        

        # ------------------------------------------------------------------------
        # Convert table into dataframe
        # ------------------------------------------------------------------------

        # converts json to dict so we can accesss the data entries 
        datadict = json.loads(beatable.text)['BEAAPI']['Results']['Data']

        # this parses the remaining dictionary entries to create a long table
        dftable = pd.DataFrame.from_records(datadict, index=range(len(datadict)))

        if Year == Yearstart:
            df = dftable
        else:
            df = df.append(dftable)


    # ------------------------------------------------------------------------
    # Format data frame to time series
    # ------------------------------------------------------------------------

    # Convert time period to datetime
    df['TimePeriod'] =  pd.to_datetime(df['TimePeriod'], format='%YM%m')

    # Set panel variables
    df.set_index(['LineDescription','TimePeriod'], inplace=True)

    # make data numeric
    df['DataValue'] = df['DataValue'].replace(',','', regex=True).astype(float)

    # convert to monthly and billions of dollars
    if TableName in ['U20405','U20406']:
        df['DataValue'] = df['DataValue'] / 12 / 1000

    # sort dataframe to allow for indexing
    df = df.sort_index(axis=0)

    # ------------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------------

    df.to_parquet('../output/' + TableTitle + '.parquet')

# %%
