## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/psmjsampleinterview.parquet \
		../input/psmjsamplemonthly.parquet \


# This generates the list of processed data files	
output=\
		../output/weightsdcartknlags0insampleinterview.parquet \
		../output/dcartknlags0insampleinterviewrbtindicatortreatment.png



# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded cex data. if statement checks that the file exists
$(input): | ../input
	if [ -e ../../psmjsample/output/$(@F) ] ; then ln -s ../../psmjsample/output/$(@F) $@ ; else exit 1; fi	

# testing for smoothing
../output/weightsdcartknlags0insampleinterview.parquet: weightsandeffects.py cohortweights.py $(input) 
	rm -f ../output/*
	python3 $<

weightfiles = $(shell find ../output/weights*.parquet)
# $(weightfiles)

../output/dcartknlags0insampleinterviewrbtindicatortreatment.eps: plotweights.py  ../output/weightsdcartknlags0insampleinterview.parquet $(weightfiles)
	python3 $< $(weightfiles)