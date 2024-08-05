## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we expect 	
output= \
		../output/ucc_category_map.yml \
		../output/level3to2_map.yml \
		../output/psmjcategories_map.yml \
		../output/cexname_conversion.yml \
		../output/nipacategories_map.yml \
		../output/ucc_nipa_map.yml \


## Recipes
# all: $(folders) $(input)  $(output)
all: $(folders) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# this is the recipe that creates the yml file mapping ucc codes to expenditure
# categories
# %: will repeat this structure for both output files
# $< is first dependency (in this case the python file)
../output/%.yml: %.py
	python $<

# ../output/level3to2_map.yml: level3to2_map.py
# 	python level3to2_map.py

