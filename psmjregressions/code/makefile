## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/psmjsampleinterview.dta  \
		../input/psmjsamplemonthly.dta  \
		../input/nipavariablesinterview.dta \
		../input/agglabelsnipa.dta \
		../input/agglabels_levelnames_nipa.dta \
		../input/agglabelspsmj.dta \


# This generates the list of output files
output= \
		../output/Table_E1_d_pce.tex \
		../output/Table_E2.tex \
		../output/Table_bias.tex \
		../output/Table_AE4.tex \
		../output/Table_AE5.tex \
		../output/Table_AE4veh.tex \
		../output/Table_AEMonthly.tex \
		../output/motorvehicle_eventstudy.pdf \
		../output/table_randomtest_d_pce.tex\
		../output/Figure_fullVrebateonly_coefstack_d_pce.pdf


# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@



../input/agglabelsnipa.dta ../input/agglabels_levelnames_nipa.dta ../input/agglabelspsmj.dta:  | ../input
	if [ -e ../../createconsumptionvariables/output/$(@F) ] ; then ln -s ../../createconsumptionvariables/output/$(@F) $@ ; else exit 1; fi	


../input/nipavariablesinterview.dta:  | ../input
	if [ -e ../../nipavariables/output/$(@F) ] ; then ln -s ../../nipavariables/output/$(@F) $@ ; else exit 1; fi	

# Create sim link to file with downloaded cex data. if statement checks that the file exists
../input/psmjsampleinterview.dta ../input/psmjsamplemonthly.dta: | ../input
	if [ -e ../../psmjsample/output/$(@F) ] ; then ln -s ../../psmjsample/output/$(@F) $@ ; else exit 1; fi	


# we first delete the output folders to make sure all data is new 
$(output):  variablelabels.do Table_E1.do Table_bias.do Table_E2.do Table_AE4.do Table_AE5.do Table_E1monthly.do Figure_fullVrebateonly.do $(input)

#2008 Tables and charts
	stata-se -b run variablelabels.do
	stata-se -b run Table_E1.do
	stata-se -b run Table_bias.do
	stata-se -b run Table_E2.do
	stata-se -b run Table_AE4.do
	stata-se -b run Table_AE5.do
	stata-se -b run Table_E1Monthly.do
	stata-se -b run Figure_fullVrebateonly.do
	stata-se -b run eventstudy_graph.do
	stata-se -b run table2_schedule_interaction.do


# rm -f ../output/*