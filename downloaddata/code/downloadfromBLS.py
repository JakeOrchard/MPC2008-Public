import yaml
import requests
import zipfile
import argparse
import os
import shutil
import argparse
from random import randint
from time import sleep

# Import yaml file directory using option passed along
parser = argparse.ArgumentParser()
parser.add_argument('--downloadopt', help='yaml file with startdate and enddate for download', required=True)
args = parser.parse_args()
yamlfile = args.downloadopt


# import YAML file with download settings
with open(yamlfile, 'r') as file:
    options_list = yaml.load(file, Loader=yaml.FullLoader)


# reading start and end dates for download from yaml file
CEX_startdate = int(options_list['CEX_startdate'])
CEX_enddate = int(options_list['CEX_enddate']) + 1

years_to_download = range(CEX_startdate, CEX_enddate)

# format that we download data in
CEX_format = options_list['CEX_format']

# extract zip file to this folder
zipextractdir = '../output/cexdownload/'

# set download agent for BLS website
# this may need to be updated, otherwise the website may block the download
session = requests.Session()
header = {
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.1 Safari/605.1.15',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
    'Accept-Encoding': 'zip, gzip',
    'Accept-Language': 'en-US,en;q=0.9',
    'Connection': 'keep-alive',
    'referer': 'https://www.bls.gov/cex/pumd_data.htm',
}

session.headers.update(header)

# this loop downloads the data from the BLS starting with CEX_startdate from 
# the yaml file and ending with CEX_enddate
for year in years_to_download:

    # last two digits of number
    yy2 = abs(year) % 100 

    # filename
    filename = "intrvw" + str(yy2).zfill(2) + ".zip"

    # construct URL (can change data format by replacing stata with comma)
    url = "https://www.bls.gov/cex/pumd/data/" + CEX_format + "/" + filename

    # construct directory
    zipfilepath = '../output/'  + filename 

    # download data
    print(url + '\n')
    # wget.download(url, zipfilepath)
    response = session.get(url)
    with open(zipfilepath, mode='wb') as localfile:     
        localfile.write(response.content)
    sleep(randint(10,100))

    # unzip all files then zip them together again. this makes it easier to track
    # the downloaded file with make
    with zipfile.ZipFile(zipfilepath, 'r') as zip_ref:
        zip_ref.extractall(zipextractdir)

        # delete downloaded zip file
        if os.name != 'nt':         #Windows throws an access error in the next line
            os.remove(zipfilepath) 

# Once all files are extracted combine all of them into a signle zip file so we 
# can track it in make
shutil.make_archive('../output/cexdownload','zip',zipextractdir)

# delete extracted files
shutil.rmtree(zipextractdir)

