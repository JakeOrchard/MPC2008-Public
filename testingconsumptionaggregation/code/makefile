## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/fmliquarterly.parquet \
		../input/mtbimonthly.parquet \
		../input/cexname_conversion.yml \
		../input/level3to2_map.yml


# This generates the list of processed data files	
output=../output/testingaggregation.txt

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to utility file
mpc_utilities.py:  | ../input
	if [ -e ../../utilities/code/$(@F) ] ; then ln -s ../../utilities/code/$(@F) $@ ; else exit 1; fi	

# Create sim link to file with downloaded cex data. if statement checks that the file exists
../input/fmliquarterly.parquet ../input/mtbimonthly.parquet: | ../input
	if [ -e ../../createconsumptionvariables/output/$(@F) ] ; then ln -s ../../createconsumptionvariables/output/$(@F) $@ ; else exit 1; fi	

# Create sim link to file with ucc mappings to expenditure categories
../input/level3to2_map.yml ../input/cexname_conversion.yml: | ../input
	if [ -e ../../ucccodemappings/output/$(@F) ] ; then ln -s ../../ucccodemappings/output/$(@F) $@ ; else exit 1; fi			

# this is the recipe that creates the zip file with cex data  --downloadopt ../input/cexdownloadoptions.yml
../output/%.txt: %.py $(input)
	rm -r ../output
	mkdir ../output
	python $<