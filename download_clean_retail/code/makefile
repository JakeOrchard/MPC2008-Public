## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# Inputs
input = \

# This generates the list of zip files we expect 	
output= ../output/retail_sales_clean.csv  \
        ../output/retail_sales_clean.dta


## Recipes
all: $(folders) $(input)  $(output)
# all: $(folders) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@


# this is the recipe that creates the zip file with cex data
../output/retail_sales_clean.csv ../output/retail_sales_clean.dta: downloadfromCensus.py 
	python downloadfromCensus.py
