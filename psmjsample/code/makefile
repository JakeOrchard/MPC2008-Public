## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/psmjvariablesinterview.parquet \
		../input/psmjvariablesmonthly.parquet 


# This generates the list of processed data files	
output= \
		../output/psmjsampleinterview.parquet \
		../output/psmjsamplemonthly.parquet \
		../output/psmjsampleinterview.dta \
		../output/psmjsamplemonthly.dta 
		


# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded cex data. if statement checks that the file exists
../input/psmjvariablesinterview.parquet ../input/psmjvariablesmonthly.parquet: | ../input
	if [ -e ../../psmjvariables/output/$(@F) ] ; then ln -s ../../psmjvariables/output/$(@F) $@ ; else exit 1; fi	

# this is the recipe that creates the zip file with cex data  --downloadopt ../input/cexdownloadoptions.yml
 ../output/psmjsampleinterview.dta ../output/psmjsampleinterview.parquet ../output/psmjsamplemonthly.dta ../output/psmjsamplemonthly.parquet: psmjsample.py $(input)
	rm -f ../output/*
	python $<
