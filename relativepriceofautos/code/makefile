## Definitions

# Folders we will be looking for in the task directory
folders := ../output  ../input

input = ../input/r-cpi-u-nv-data.xlsx

# This generates the list of processed data files and figures
output=  ../output/newvehiclerelpresearch.eps
		 

# Recipes
all: $(folders) $(output) $(input)

# Create folders if they do not already exist
$(folders):
	mkdir $@

$(input): | ../input
	if [ -e ../../external_data/$(@F) ] ; then ln -s ../../external_data/$(@F) $@ ; else exit 1; fi	

$(output): relativepriceofautos.do $(input)
	stata-se -b run $<
