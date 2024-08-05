import pandas as pd
import numpy as np
from linearmodels import PanelOLS

def cohortweights(dfin, outcome, treatment, relativetime, cohort,
                  controls=[], entity_effects=True, time_effects=True,
                  weights=None):

    # new dataframe with the variables we need
    # this way we do not override the existing dataframe
    df = dfin[[outcome, treatment, relativetime, cohort] + controls].copy()

    # add constant to datframe
    df['__cons'] = 1

    # partial out controls wrt fixed effects
    controlspartial = []
    for variable in controls:

        var_reg = PanelOLS(dependent=df[variable], 
                           exog=df[['__cons']], 
                           entity_effects=entity_effects, 
                           time_effects=time_effects,
                           weights=weights,
                           ).fit()

        df['__' + variable + 'partial'] = var_reg.resids

        controlspartial.append('__' + variable + 'partial')

    # get double-demeaned treatment and outcome variable partialled out wrt controls
    for variable in outcome, treatment:

        var_reg = PanelOLS(dependent=df[variable], 
                           exog=df[['__cons'] + controlspartial], 
                           entity_effects=entity_effects, 
                           time_effects=time_effects,
                           weights=weights,
                           ).fit()

        df['__' + variable + 'partial'] = var_reg.resids
    
    
    # Two-way FE regressions to check weighted effect calculation
    direct_est = PanelOLS(dependent=df['__' + outcome + 'partial'], 
                          exog=df[['__' + treatment + 'partial']],
                          weights=weights,
                          ).fit().params[0]                         
    
    # get groups of relative time and cohort
    time = df.index.names[1]
    dfgroupby = df.groupby([relativetime,cohort,time])
    df['__Groups'] = dfgroupby.ngroup()

    dfgroupbytest = df.groupby([time])
    df['__GroupsTime'] = dfgroupbytest.ngroup()

    # create dummies for groups
    group_dummies = pd.get_dummies(df['__Groups'], prefix='__Group_Num')

    # mapping group numbers to group size, cohort and relative time
    dfgroup = dfgroupby['__Groups'].first().to_frame()
    dfgroup = dfgroup.join(dfgroupby.size().rename('GroupSize'))
    dfgroup = dfgroup.reset_index().set_index('__Groups')
    
    # Regression of groups on treatment variable to get weights
    dfgroup['Weights'] = 0

    for g in range(len(dfgroup)):

        # do not need controls here because RHS is partialled out
        bin_reg = PanelOLS(dependent=group_dummies['__Group_Num_' + str(g)], 
                           exog=df[['__' + treatment + 'partial']],
                           weights=weights,
                           ).fit()    

        dfgroup.loc[g, 'Weights'] = bin_reg.params[0]

    # treatment effect calculated from regression partialled y on groups and controls
    y_reg = PanelOLS(dependent=df['__' + outcome + 'partial'], 
                     exog=group_dummies.join(df[controlspartial]),
                     weights=weights,
                     ).fit()  

    dfgroup['Cohort Treatment'] = 0
    
    for g in range(len(dfgroup)):
        dfgroup.loc[g, 'Cohort Treatment'] = y_reg.params['__Group_Num_' + str(g)]

    # calculate weighted average
    dfgroup['Contribution'] = dfgroup['Weights'] * dfgroup['Cohort Treatment']

    weight_est = dfgroup['Contribution'].sum()

    if abs(direct_est-weight_est)>10 ** -10:
        print('Warning: estimates do not add up')
        print(direct_est, weight_est)  

    return weight_est, dfgroup


def cohortweightsmult(dfin, outcome=[], treatment=[], relativetime=[], cohort=[],
                  controls=[], entity_effects=True, time_effects=True,
                  weights=None):   

    allvars = outcome + treatment + relativetime + cohort + controls

    # new dataframe with the variables we need
    # this way we do not override the existing dataframe
    df = dfin[allvars].copy()

    # add constant to datframe
    df['__cons'] = 1

    # Weights
    if weights!=None:
        weights = dfin[weights].copy()

    # partial out controls wrt fixed effects
    controlspartial = []
    for variable in controls:

        var_reg = PanelOLS(dependent=df[variable], 
                           exog=df[['__cons']], 
                           entity_effects=entity_effects, 
                           time_effects=time_effects,
                           weights=weights,
                           ).fit()

        df['__' + variable + 'partial'] = var_reg.resids

        controlspartial.append('__' + variable + 'partial')

    # get double-demeaned treatment and outcome variable partialled out wrt controls
    for variable in outcome + treatment:

        var_reg = PanelOLS(dependent=df[variable], 
                           exog=df[['__cons'] + controlspartial], 
                           entity_effects=entity_effects, 
                           time_effects=time_effects,
                           weights=weights,
                           ).fit()

        df['__' + variable + 'partial'] = var_reg.resids
    
    
    # Two-way FE regressions to check weighted effect calculation
    yvar = '__' + outcome[0] + 'partial'
    treat_list = ['__' + treat + 'partial' for treat in treatment]

    direct_est = PanelOLS(dependent=df[yvar], 
                          exog=df[treat_list],
                          weights=weights,
                          ).fit().params[treat_list]                         

    # get groups of relative time and cohort
    time = [df.index.names[1]]
    dfgroupby = df.groupby(relativetime + cohort + time)
    df['__Groups'] = dfgroupby.ngroup()

    dfgroupbytest = df.groupby(time)
    df['__GroupsTime'] = dfgroupbytest.ngroup()

    # create dummies for groups
    group_dummies = pd.get_dummies(df['__Groups'], prefix='__Group_Num')

    # mapping group numbers to group size, cohort and relative time
    dfgroup = dfgroupby['__Groups'].first().to_frame()
    dfgroup = dfgroup.join(dfgroupby.size().rename('GroupSize'))
    dfgroup = dfgroup.reset_index().set_index(['__Groups'] + relativetime + cohort + time + ['GroupSize'])

    # creating dataframe for estimation results and merging group info
    group_estimation_index = pd.MultiIndex.from_product(
                                 [list(range(len(group_dummies.columns))),treatment],
                                 names = ['__Groups', 'Variable']
                                 )

    dfestimates = pd.DataFrame(index = group_estimation_index,
                               columns = ['Weights', 'Cohort Treatment'],
                               dtype = float)

    dfestimates = dfestimates.merge(dfgroup, how='left', left_index=True, right_index=True)                               


    for g in range(len(dfgroup)):

        # do not need controls here because RHS is partialled out
        bin_reg = PanelOLS(dependent=group_dummies['__Group_Num_' + str(g)], 
                           exog=df[treat_list],
                           weights=weights,
                           ).fit()    

        for treat_var, treat_name in zip(treat_list, treatment):
            idx = pd.IndexSlice[g, treat_name]
            dfestimates.loc[idx, 'Weights'] = bin_reg.params[treat_var]

    # treatment effect calculated from regression partialled y on groups and controls
    for treat_var, treat_name in zip(treat_list, treatment):
        
        y_reg = PanelOLS(dependent=df[yvar], #yvar -- check individual levels
                        exog=group_dummies.join(df[controlspartial]),
                        weights=weights,
                        ).fit()  


        for g in range(len(dfgroup)):
            idx = pd.IndexSlice[g, treat_name]
            dfestimates.loc[idx, 'Cohort Treatment'] = y_reg.params['__Group_Num_' + str(g)]  


    dfestimates['Contribution'] = dfestimates['Weights']*dfestimates['Cohort Treatment']  

    # Check aggregation
    check_aggregation(dfestimates, direct_est, treatment)           

    # merge group characteristics
    dfestimates = dfestimates.merge(dfgroup, how='left', left_index=True, right_index=True)

    return dfestimates    


def check_aggregation(dfestimates, direct_est, treatment):
    
    dfsum = dfestimates['Contribution'].groupby('Variable').sum()
    
    for var in treatment:
        assert np.abs(dfsum[var] - direct_est['__' + var + 'partial'])<10**-8
