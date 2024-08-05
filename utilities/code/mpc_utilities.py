def aggregate_with_dict(df,filename):  
    # function aggregates df using the dictionary contained in filename
    # assumes that the file is stored in yaml format

    import yaml 

    # import YAML file with dictionary that maps several categories into one
    with open(filename, 'r') as yamlfile:
        dictionary_map = yaml.load(yamlfile, Loader=yaml.FullLoader)

    # loop over categories and sum entries
    for category, subcategories in dictionary_map.items():
        df[category] = df[subcategories].sum(axis=1)

    return df


def unzip_import_csv_as_df(filename, cols_w_index=None, cols_to_import=None):

    print('USE ONLY FOR FILES TO BE IMPORTED INTO STATA. OTHERWISE USE PARQUET')

    import zipfile
    import shutil
    import pandas as pd

    # directory we import from:
    zipimportdir = '../input/'

    # directory we extract to:
    zipextractdir = '../output/' + filename + '/'
    
    # unzip files
    with zipfile.ZipFile(zipimportdir + filename + '.zip', 'r') as zip_ref:
            zip_ref.extractall(zipextractdir)

    # load csv file as dataframe
    df = pd.read_csv(zipextractdir + filename + '.csv',
                        index_col   = cols_w_index,
                        usecols     = cols_to_import) 

    # delete directory with extracted files    
    shutil.rmtree(zipextractdir)  

    return df    

def save_df_csv_and_zip(df, filename, save_index=True):

    print('USE ONLY FOR FILES TO BE IMPORTED INTO STATA. OTHERWISE USE PARQUET')

    import zipfile
    import os
    import pandas as pd

    # directory we extract to:
    zipsavename = '../output/' + filename 

    # save combined file as csv
    df.to_csv(zipsavename + '.csv', index=save_index)
    
    # and zip the file
    with zipfile.ZipFile(zipsavename + '.zip', 'w', zipfile.ZIP_DEFLATED) as zip_write:
            zip_write.write(zipsavename + '.csv', arcname=filename +'.csv')

    # delete csv file
    os.remove(zipsavename + '.csv')  

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