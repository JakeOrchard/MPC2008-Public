## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input = \
		../input/pceexpenditure.parquet \
		../input/pceexpenditureprice.parquet \
		../input/pceexpenditurereal.parquet \
		../input/forecasts.dta \
		../input/forecasts_lowerCI.dta \
		../input/freddata_for_forecasting.dta \
		../input/motorvehicleprice.dta
	

inputrebate = ../input/rel_rebate.dta \

modelfiles = nk_rebates_mpc_model.py \
			 nk_rebates_mpc_monthly_ss.py

# This generates the list of processed data files	
output= ../output/mpcsgammadurableprice.tex 
outputreg = ../output/modelregression.tex 

# Recipes
all: $(folders) $(input) $(inputjps) $(inputrebate) $(output) $(outputreg) 
# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded cex data. if statement checks that the file exists
../input/pceexpenditure.parquet ../input/pceexpenditurereal.parquet ../input/pceexpenditureprice.parquet: | ../input
	if [ -e ../../downloaddata/output/$(@F) ] ; then ln -s ../../downloaddata/output/$(@F) $@ ; else exit 1; fi	

../input/forecasts.dta ../input/forecasts_lowerCI.dta ../input/freddata_for_forecasting.dta: | ../input
	if [ -e ../../forecasting/output/$(@F) ] ; then ln -s ../../forecasting/output/$(@F) $@ ; else exit 1; fi		

../input/motorvehicleprice.dta: | ../input
	if [ -e ../../relativepriceofautos/output/$(@F) ] ; then ln -s ../../relativepriceofautos/output/$(@F) $@ ; else exit 1; fi		


$(inputrebate): | ../input
	if [ -e ../../narrative/output/$(@F) ] ; then ln -s ../../narrative/output/$(@F) $@ ; else exit 1; fi		

# run model file and move outputs
$(output): nk_rebates_mpc_analysis.py nk_rebates_mpc_plot.py $(input) $(modelfiles) $(inputrebate)
	rm -f ../output/*
	python $<
	python nk_rebates_mpc_plot.py

# run model regressions
 $(outputreg): nk_rebates_mpc_regression.py modelregressions.do $(modelfiles)
	python $<	
	stata-se -b run modelregressions.do

# baseline_irfs.xlsx: rot_loopm.m nk_rebates_mpc_monthly_loops.mod  $(input) 
# 	rm -f ../output/*
# 	matlab -nojvm -nodisplay -nosplash < "rot_loopm.m"
# 	mv nk_rebates_mpc_monthly_loops.log paramfile.mat nk_rebates_mpc_monthly_loops_results.mat ../output
# 	rm -r nk_rebates_mpc_monthly_loops +nk_rebates_mpc_monthly_loops
	

# ../output/mpc_counterfactual.eps: rebate2008_counterfactual.do baseline_irfs.xlsx
# 	stata-se -b run rebate2008_counterfactual.do
# 	mv rebate2008_counterfactual.log ../output


# ../output/mpc_counterfactual2.eps: nk_rebates_mpc_analysis.py ../input/pceexpenditure.parquet
# 	python nk_rebates_mpc_analysis.py
