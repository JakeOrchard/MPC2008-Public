## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/pceexpenditure.parquet


# This generates the list of processed data files	
output=\
		../output/counterfactuals.eps



# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded cex data. if statement checks that the file exists
$(input): | ../input
	if [ -e ../../downloaddata/output/$(@F) ] ; then ln -s ../../downloaddata/output/$(@F) $@ ; else exit 1; fi	

# testing for smoothing
../output/counterfactuals.eps: pcestatistics.py $(input) 
	rm -f ../output/*
	python $<