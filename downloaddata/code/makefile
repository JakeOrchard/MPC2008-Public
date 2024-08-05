## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# Inputs
input = ../input/cexdownloadoptions.yml

# This generates the list of zip files we expect 	
output=../output/cexdownload.zip  \
		../output/pceexpenditure.parquet


## Recipes
all: $(folders) $(input)  $(output)
# all: $(folders) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to options file. if statement checks that the file exists
../input/cexdownloadoptions.yml: | ../input
	if [ -e ../../globaloptions/hand/$(@F) ] ; then ln -s ../../globaloptions/hand/$(@F) $@ ; else exit 1; fi	

# this is the recipe that creates the zip file with cex data
../output/cexdownload.zip: downloadfromBLS.py ../input/cexdownloadoptions.yml
	python downloadfromBLS.py --downloadopt ../input/cexdownloadoptions.yml

# this is the recipe that creates the zip file with cex data  --downloadopt ../input/cexdownloadoptions.yml
../output/pceexpenditure.parquet: pcefromBEA.py
	python pcefromBEA.py


# This is the link to ICPSR but it requires log-in and agreement. R has a package called icpsrdata that could automate this process
# http://www.icpsr.umich.edu/cgi-bin/bob/zipcart2?path=ICPSR&study=8423&bundle=all&ds=&dups=yes