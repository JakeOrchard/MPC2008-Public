#%%
import pickle
from re import X
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from dateutil.relativedelta import relativedelta
from cycler import cycler


# ------------------------------------------------------------------------
# Tables and Figure Settings
# ------------------------------------------------------------------------

parasetnames = ['baseline'] 

info_set = ['no-anticipation']

table_names = ['PCE', 'Motor vehicles', 'Nondurable goods'] #,

model_names = dict()
model_names['Real'] = ['exp', 'realx', 'realc']
model_names['Nominal'] = ['nomexp', 'nomx', 'nomc']

# graph defaults
plt.rc('font', size=12)
new_prop_cycle = cycler('color', ['k','tab:blue','tab:orange','tab:green','r'])
new_prop_cycle_model = cycler('color', ['tab:blue','tab:orange','tab:green','r'])
plt.rc('axes', prop_cycle=new_prop_cycle)

# ------------------------------------------------------------------------
# Regression settings
# ------------------------------------------------------------------------

# interview sample size in each month
n = 4500
rbtamt = 950

# fraction treated each month, 0 = April rebate
# 0,3,6 on same cycle
treat = np.zeros([9,])
treat[0] = 0.030
treat[1] = 0.266
treat[2] = 0.279
treat[3] = 0.141
treat[4] = 0.026
treat[5] = 0.011
treat[6] = 0.006
treat[7] = 0.005
treat[8] = 0.002
# treat[9] = 0.003
# treat[10] = 0
# treat[11] = 0

# length of IRF
hhlen = 24

# ------------------------------------------------------------------------
# Tables and Figures
# ------------------------------------------------------------------------


for parasetname in parasetnames:

    output = pickle.load( open( '../output/' + parasetname + '.pkl', 'rb' ) )
    parameter_set = output.keys()

    for parameter in parameter_set:

        mpcset = output[parameter].keys()

        for mpc in mpcset:

            ss = output[parameter][mpc]['steady']

            G = output[parameter][mpc]['G']['micro']

            Tirf, Tnews = G['cr']['incshock'].shape

            # nondurable exp
            for hh in ['r','o']:
                G['nd' + hh] = dict()
                G['nd' + hh]['incshock'] = G['c' + hh]['incshock'] + ss['eta']*G['d' + hh]['incshock']

            # 3 - month sum of lagged expenditures
            for exp in ['xr','ndr','xo','ndo']:
                G['3' + exp] = dict()
                G['3' + exp]['incshock'] = (  G[exp]['incshock'] 
                                        +   np.vstack([np.zeros([1,Tnews]),G[exp]['incshock'][:-1,:]])
                                        +   np.vstack([np.zeros([2,Tnews]),G[exp]['incshock'][:-2,:]])
                                        )
            
            


            df = pd.DataFrame()

            

            dates = np.arange(0, hhlen)
            ndf = 0

            for treattime, probtreat in enumerate(treat):

                for hh, frachh in zip(['r','o'],[np.round(ss['gamma'],2),np.round(1-ss['gamma'],2)]):
                    data = np.zeros([4 + hhlen,6])
                    for i,var in enumerate(['3x','3nd']):
                        levelc = np.zeros([7 + hhlen,])
                        levelc[7:] = G[var + hh]['incshock'][:hhlen,2+treattime]  #-  G[var + hh]['incshock'][0: hhlen,3+treattime]
                        data[:,i*2] = levelc[3:]
                        # diffc = levelc[3:] - levelc[:-3]
                        data[:,i*2+1] = levelc[3:] - levelc[:-3]
                        
                    data = data * rbtamt
                    
                    
                    
                    for inttime in range(3):

                        inttimeshift = 2

                        for treatintnum in range(3,6):
                           

                            datanum = data.copy()

                            datanum = datanum[inttimeshift + treattime + (5-treatintnum)*3:inttimeshift + treattime + (5 - treatintnum)*3 + 9:3]

                            datanum[treatintnum-3,4] = 1
                            datanum[treatintnum-3,5] = rbtamt


                            Tpanel, trash = datanum.shape

                            # for datatreat, fracpop in zip([datanum,datacontrol], [probtreat,1-probtreat]):
                            nhh = np.int32(np.round(n * frachh * probtreat))
                            idhh = np.tile(ndf + np.arange(0, nhh), (Tpanel,1)).reshape(nhh*Tpanel, 1, order='F')
                            datevec = np.tile(dates[6+treattime+inttime + (3-treatintnum)*3:6+treattime+inttime + (3 - treatintnum)*3 + 9:3,np.newaxis], (nhh,1))
                            numvec = np.tile(np.array([3,4,5])[:,np.newaxis], (nhh,1))

                            dataset = np.tile(datanum, (nhh,1))
                            index = pd.MultiIndex.from_arrays(np.concatenate((idhh,datevec,numvec), axis=1).T, names=('Household','Date','IntNum'))

                            dftemp = pd.DataFrame(data=dataset, index=index, columns=['Durable','d_Durable','Nondurable','d_Nondurable','Rebate','RebateAmount'])
                            df = df.append(dftemp)
                            ndf = ndf + nhh

                        
                            
            
            

            for date in range(df.index.get_level_values('Date').min(),df.index.get_level_values('Date').max() + 1):
                
                if df.loc[pd.IndexSlice[:,date,:],'Rebate'].sum()==0:
                    continue

                for intnum in range(3,6):
                    
                    # count how many control households we will need 
                    count = df.loc[pd.IndexSlice[:,date,intnum],'Rebate'].count()
                    print(date, count)

                    # control group
                    datanum = np.zeros([3,6])

                    nhh = np.int32(n - count)
                    idhh = np.tile(ndf + np.arange(0, nhh), (3,1)).reshape(nhh*3, 1, order='F')

                    datevec = np.tile(dates[date + (3-intnum)*3:date + (3 - intnum)*3 + 9:3,np.newaxis], (nhh,1))
                    numvec = np.tile(np.array([3,4,5])[:,np.newaxis], (nhh,1))

                    dataset = np.tile(datanum, (nhh,1))
                    index = pd.MultiIndex.from_arrays(np.concatenate((idhh,datevec,numvec), axis=1).T, names=('Household','Date','IntNum'))

                    dftemp = pd.DataFrame(data=dataset, index=index, columns=['Durable','d_Durable','Nondurable','d_Nondurable','Rebate','RebateAmount'])
                    
                    df = df.append(dftemp)

                    ndf = ndf + nhh
                    

            
            df.to_stata('../output/simul' + parasetname + parameter + 'mpc' + mpc[-3:].replace('.','') + '.dta')
# %%
