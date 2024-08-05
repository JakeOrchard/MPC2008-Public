## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/nipavariablesmonthly.dta \
		../input/retail_sales_clean.dta \

# This generates the list of output figures
cexsa = ../output/sa_cex.dta

output= \
		../output/furn_cexvnipa.pdf \
		../output/durables_cexvnipa.pdf \
		../output/mvparts_cexvnipa.pdf \
		../output/nondurables_cexvnipa.pdf \
		../output/food_cexvnipa.pdf \
		../output/clothing_cexvnipa.pdf \
		../output/gas_cexvnipa.pdf \
		../output/services_cexvnipa.pdf \
		../output/pce_cexvnipa.pdf \
		../output/pce_cex_series.dta \
		../output/billion_difference.pdf \
		../output/log_difference.pdf \
		../output/fig_real_rts_pcegoods.pdf 

		
# Recipes
all: $(folders) $(input) $(output) $(cexsa)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to input files 
../input/nipavariablesmonthly.dta: | ../input
	if [ -e ../../nipavariables/output/$(@F) ] ; then ln -s ../../nipavariables/output/$(@F) $@ ; else exit 1; fi	

../input/retail_sales_clean.dta: | ../input
	if [ -e ../../download_clean_retail/output/$(@F) ] ; then ln -s ../../download_clean_retail/output/$(@F) $@ ; else exit 1; fi

# Recipe for output files
#../output/sa_cex.dta: season_adjust_cex.do ../input/nipavariablesmonthly.dta 
#	stata-se -b run season_adjust_cex.do

#../output/furn_cexvnipa.pdf ../output/durables_cexvnipa.pdf ../output/mvparts_cexvnipa.pdf 	../output/nondurables_cexvnipa.pdf 	../output/food_cexvnipa.pdf ../output/clothing_cexvnipa.pdf ../output/gas_cexvnipa.pdf ../output/services_cexvnipa.pdf ../output/pce_cexvnipa.pdf ../output/pce_cex_series.dta: graph_pceVcex.do ../output/sa_cex.dta 
#	stata-se -b run graph_pceVcex.do


#../output/billion_difference.pdf ../output/log_difference.pdf : dist_measurement_error.do ../output/pce_cex_series.dta 
#	stata-se -b run dist_measurement_error.do


$(cexsa): season_adjust_cex.do  $(input)
	rm -f ../output/*
	stata-se -b run  season_adjust_cex.do

$(output): graph_pceVcex.do graph_retailVcexVPCE.do dist_measurement_error.do   $(input) $(cexsa)
	stata-se -b run graph_pceVcex.do
	stata-se -b run graph_retailVcexVPCE.do
	stata-se -b run dist_measurement_error.do	

# rm *.log
