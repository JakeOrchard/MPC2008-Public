## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of processed data files and figures
input = ../input/sca-tableall-on-2021-Dec-02.csv

output=  ../output/clean_michigan.dta\
		 

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded michigan
../input/sca-tableall-on-2021-Dec-02.csv: | ../input
	if [ -e ../../external_data/$(@F) ] ; then ln -s ../../external_data/$(@F) $@ ; else exit 1; fi	


# this is the recipe that creates the clean michigan data
../output/clean_michigan.dta: clean_michigan.do $(input)
	stata-se -b run clean_michigan.do
