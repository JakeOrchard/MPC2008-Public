## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of processed data files and figures
input = ../input/levels_SPF0704.dta\
		../input/greenbook08.dta \
		../input/forecasts.dta  

output=  ../output/all_forecasts.eps 
		 

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Link Inputs
../input/greenbook08.dta: | ../input
	if [ -e ../../greenbook_forecast/output/$(@F) ] ; then ln -s ../../greenbook_forecast/output/$(@F) $@ ; else exit 1; fi	

../input/levels_SPF0704.dta: | ../input
	if [ -e ../../survey_prof_fore/output/$(@F) ] ; then ln -s ../../survey_prof_fore/output/$(@F) $@ ; else exit 1; fi	

../input/forecasts.dta: | ../input
	if [ -e ../../forecasting/output/$(@F) ] ; then ln -s ../../forecasting/output/$(@F) $@ ; else exit 1; fi	


# this is the recipe 
$(output): graph_forecasts_all.do $(input)
	stata-se -b run graph_forecasts_all.do
