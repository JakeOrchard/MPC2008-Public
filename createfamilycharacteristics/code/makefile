## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/fmliquarterly.parquet 


# This generates the list of processed data files	
output=\
		../output/cohortvariables.parquet \
		../output/interviewvariables.parquet


# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded cex data. if statement checks that the file exists
../input/fmliquarterly.parquet: | ../input
	if [ -e ../../createconsumptionvariables/output/$(@F) ] ; then ln -s ../../createconsumptionvariables/output/$(@F) $@ ; else exit 1; fi	

# this is the recipe that creates the zip file with cex data  --downloadopt ../input/cexdownloadoptions.yml
 ../output/%.parquet: %.py $(input)
	python $<
