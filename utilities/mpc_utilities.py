import re
import yaml 
import pandas as pd

def aggregate_with_dict(df,filename):  
    """
    function aggregates df using the dictionary contained in filename
    assumes that the file is stored in yaml format
    """
    

    # import YAML file with dictionary that maps several categories into one
    with open(filename, 'r') as yamlfile:
        dictionary_map = yaml.load(yamlfile, Loader=yaml.FullLoader)

    # loop over categories and sum entries
    for category, subcategories in dictionary_map.items():
        df[category] = df[subcategories].sum(axis=1)

    return df


def aggregate_df(df, agg_by={}, **kargs):
    """
    This function takes a dataframe, and aggregates variables using the operation specified in positional arguments.
    For example, to sum a list of varibales write (also accepts set/tuple input):
    sum = [list of variables]
    
    The variables in agg_by will be in the index of the returned dataframe. The variables in the columns are those 
    specified in kargs.
    """    
    # return original dataframe if empty set
    if not agg_by:
        print('aggregate_df: No Aggregation Variables. Returning original DF.')
        return df
    
    # first we map each variable to the operation and collect the variables
    # e.g. sum=var1 gets mapped to 'var1':'sum'
    aggregation_dict = {var:operation for operation, variables in kargs.items() for var in variables} 
        
    # new dataframe: combines all series in agg_by and in kargs
    df = df[set(agg_by) | aggregation_dict.keys()]
    
    # aggregate
    df = df.groupby(list(agg_by)).agg(aggregation_dict)
    
    return df    


def savetexscalar(value, name, decimals=2):
    """
    Saves value to tex file so it can be called from a latex document with
    \[name]math in math mode and \[name]text in text mode.
    Decimals are the maximum number of decimals printed. Trailing 0s and decimal
    points will not be printed.

    Note: names cannot contain numbers or underscores for latex to work
    """

    name_string = str(name)
    if re.search(r'[0-9]|\_+',  name_string):
        print('Warning: numbers or underscores not supported as Latex locals.')

    format_string = '%.' + str(decimals) + 'f'
    formatted_value = str((format_string % value).rstrip('0').rstrip('.'))
    f = open('../output/' + name_string + '.tex', 'w')
    f.write('\\newcommand{\\' + name_string + 'math}{' + formatted_value + '} \n')
    f.write('\\newcommand{\\' + name_string + 'text}{\\textnormal{' + formatted_value + ' }} \n') 
    f.close()


def fill_missing_by_type(df):
    """
    Fills missing values in a dataframe by type:
    0s for floats and integers
    False for boolean and objects
    2200-01-01 for date variables
    """
    for var in df.columns:

        if df[var].dtype in ['int64', 'float64']:
            df[var] = df[var].fillna(0)

        elif df[var].dtype in ['bool', 'object']:     
            df[var] = df[var].fillna(False) 

        elif df[var].dtype in ['datetime64[ns]']: 
            df[var] = df[var].fillna(pd.to_datetime('2200-01-01')) 

    return df   
