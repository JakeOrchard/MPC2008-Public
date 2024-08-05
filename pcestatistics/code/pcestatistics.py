#%%
import pandas as pd
import datetime
import pandas_datareader.data as web  
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------
# Load recession indicators and PCE data
# ------------------------------------------------------------------------



start = datetime.datetime(1959, 1, 1)
end = datetime.datetime(2020, 12, 31)
recession = web.DataReader(["USRECD"], "fred", start, end)

#%%
df = pd.read_parquet('../input/pceexpenditure.parquet')

# ------------------------------------------------------------------------
# Dataseries for which to compute summary stats
# ------------------------------------------------------------------------

expenditure_series = ['Personal consumption expenditures', 'New motor vehicles (55)', 'Net purchases of used motor vehicles (56)']

# create time series dataframe
df = df['DataValue']
df = df.loc[expenditure_series]
df = df.reorder_levels(order=['TimePeriod','LineDescription']).unstack()
df = df.merge(recession, how='left', left_index=True, right_index=True)


# ------------------------------------------------------------------------
# Create variables
# ------------------------------------------------------------------------

dfgrowth = pd.DataFrame(index=df.index)
dfgrowth['USRECD'] = df['USRECD']
for var in expenditure_series:
    for diff in range(4):
        dfgrowth[str(diff) + '-Month Growth ' + var] = (df[var] / df[var].shift(periods=diff, freq='MS') - 1)*100

#%%
# ------------------------------------------------------------------------
# Implied consumption paths
# ------------------------------------------------------------------------

# print(dfgrowth.describe())

# print(dfgrowth.quantile(0))
# print(dfgrowth.quantile(0.05))

growthlist = [col for col in dfgrowth.columns if col!='USRECD']

enddatelist = ['2008-08-31','2019-12-31']

for var in expenditure_series:
    dfcounterfactual = pd.DataFrame(index=df.index)
    dfcounterfactual[var] = df[var]
    dfcounterfactual = dfcounterfactual['2008-03-01':'2008-08-31']    

    results = dict()

    for enddate in enddatelist:
        print(enddate)
        dfgrowthrec = dfgrowth[:enddate]
        dfgrowthrec = dfgrowthrec[dfgrowthrec['USRECD']==1]
        # print(dfgrowthrec.describe())

        # print('Consumption Expenditure Growth Conditional on Recession')
        for quantile in [0, 0.05, 0.1, 0.25, 0.5]:
            # print('Percentile: ' + str(quantile*100) + '%')
            # print(dfgrowth[dfgrowth['USRECD']==1][growthlist].quantile(quantile))
            # print('\n')

            
            varcount = 'Percentile ' + str(quantile*100)

            for months in range(4):
                newdate = pd.to_datetime('2008-05-01') + pd.DateOffset(months=months)
                growthvar = str(months) + '-Month Growth ' + var

                dfcounterfactual.loc[newdate, varcount] = dfcounterfactual.loc['2008-05-01', var] * (1 + dfgrowthrec[growthvar].quantile(quantile)/100)


        results[enddate] = dfcounterfactual.copy()


    plt.figure()
    ax1 = plt.subplot(2,1,1)
    h1 = results[enddatelist[0]].plot(ax=ax1, legend=None)
    ax1.set_title(var + ' data through ' + enddatelist[0])
    plt.legend(bbox_to_anchor=(1,1), loc="upper left")

    ax2 = plt.subplot(2,1,2, sharex=ax1)
    h2 = results[enddatelist[1]].plot(ax=ax2, legend=None)
    ax2.set_title(var + ' data through ' + enddatelist[1])
    
    
    print('pcestatistics: update by category')
    plt.savefig('../output/counterfactuals.eps', bbox_inches='tight')
# plt.legend( handles=[h1, h2],loc="upper left", bbox_to_anchor=[0, 1])
# ncol=2, shadow=True, title="Legend", fancybox=True)


# %%
