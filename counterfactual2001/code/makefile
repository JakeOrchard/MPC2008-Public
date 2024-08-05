## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/forecasts2001.dta  \
		../input/forecasts2001_errorbands.dta \
		../input/JPS_consumption_rebate.xlsx \



# This generates the list of output files
output= \
		../output/baseline_2001_irfs.csv \
		../output/fig_jpscat_cf.eps \
		../output/fig_jpscat_gecf.eps \
		../output/fig_jpscat_gecf_for.eps \

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

../input/JPS_consumption_rebate.xlsx: | ../input
	if [ -e ../../external_data/$(@F) ] ; then ln -s ../../external_data/$(@F) $@ ; else exit 1; fi	

../input/forecasts2001.dta ../input/forecasts2001_errorbands.dta:  | ../input
	if [ -e ../../forecasting/output/$(@F) ] ; then ln -s ../../forecasting/output/$(@F) $@ ; else exit 1; fi	

# run model file and move outputs
$(output):  calldynare.m jps_counterfactual.do $(input)
	rm -f ../output/*
	matlab -r calldynare
	stata-se -b run jps_counterfactual.do