## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output 

# This generates the list of processed data files and figures

input= \
	 ../input/Motor_vehicle_data.xlsx \
	 ../input/rebates.xlsx \
	 ../input/retail_pce_dat.xlsx \
	 ../input/Update_of_Sahm-Shapiro-Slemrod-Motorvehicles.xlsx \

outputdata =  	   ../output/narrative_data.dta

outputcounter =  ../output/fig_sss_mv_counter.eps \
				 ../output/fig_sss_mv_counter_alt.eps 

outputnarrative =  ../output/fig_cy_combo.eps \
				   ../output/fig_rebates.eps \
				   ../output/fig_pce_type.eps \
				   ../output/fig_pconsnorm_combo.eps \
				   ../output/fig_ffrante_combo.eps 

outputmotor = ../output/fig_mv_units_combo.eps		 

# Recipes
all: $(folders) $(input) $(outputdata) $(outputnarrative) $(outputcounter) $(outputmotor)


# Create folders if they do not already exist
$(folders):
	mkdir $@

../input/Motor_vehicle_data.xlsx	 ../input/Update_of_Sahm-Shapiro-Slemrod-Motorvehicles.xlsx ../input/retail_pce_dat.xlsx ../input/rebates.xlsx :  | ../input
		if [ -e ../../external_data/$(@F) ] ; then ln -s ../../external_data/$(@F) $@ ; else exit 1; fi	

../input/pce_apparel_health.csv:  | ../input
		if [ -e ../../download_NIPA/output/$(@F) ] ; then ln -s ../../download_NIPA/output/$(@F) $@ ; else exit 1; fi	


$(outputdata): build_narrative_data.do $(input)
	stata-se -b run $<

$(outputnarrative): narrative_graphs.do $(outputdata)
	stata-se -b run $<	

$(outputcounter): psjm_sss_counter.do ../input/Update_of_Sahm-Shapiro-Slemrod-Motorvehicles.xlsx | ../input 
	stata-se -b run $<	

$(outputmotor): motor_vehicle_graphs.do ../input/Motor_vehicle_data.xlsx | ../input 
	stata-se -b run $<	