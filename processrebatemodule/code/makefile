## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/rbt.parquet \
		../input/fmliquarterly.parquet \
		../input/mtbimonthly.parquet \
		../input/TAX.parquet 

# This generates the list of processed data files	
output= \
		../output/rebatecuid.parquet \
		../output/rebatemonthly.parquet \
		../output/rebateinterview.parquet \




# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file RBT file
../input/rbt.parquet ../input/TAX.parquet: | ../input
	if [ -e ../../appendCEXfiles/output/$(@F) ] ; then ln -s ../../appendCEXfiles/output/$(@F) $@ ; else exit 1; fi	


# Create sim link to file FMLI file
../input/fmliquarterly.parquet ../input/mtbimonthly.parquet: | ../input
	if [ -e ../../createconsumptionvariables/output/$(@F) ] ; then ln -s ../../createconsumptionvariables/output/$(@F) $@ ; else exit 1; fi	

# this is the recipe that creates the rebate files
../output/rebatecuid.parquet ../output/rebatemonthly.parquet ../output/rebateinterview.parquet: rebate.py $(input)
	rm -f ../output/*
	python $<


