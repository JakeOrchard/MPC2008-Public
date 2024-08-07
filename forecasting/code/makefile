## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= ../input/clean_michigan_small.dta \
		../input/auxiliary_forecasting_data.csv\
		../input/rebates.csv \
		../input/JPS_consumption_rebate.xlsx


# This generates the list of processed data files	
output= ../output/forecasts.dta  \


# Recipes
all: $(folders) $(output) $(input)


# Create folders if they do not already exist
$(folders):
	mkdir $@

../input/clean_michigan_small.dta: | ../input
	if [ -e ../../michigan_survey/output/$(@F) ] ; then ln -s ../../michigan_survey/output/$(@F) $@ ; else exit 1; fi	

../input/auxiliary_forecasting_data.csv	../input/rebates.csv ../input/JPS_consumption_rebate.xlsx: | ../input
	if [ -e ../../external_data/$(@F) ] ; then ln -s ../../external_data/$(@F) $@ ; else exit 1; fi	


# run model file and move outputs
$(output): forecast.do build_forecast_data.do $(input)
	stata-se -b run build_forecast_data.do
	stata-se -b run $<

	

