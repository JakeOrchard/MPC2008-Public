## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of processed data files and figures
input = 

output= ../output/spf_rcons.dta\
../output/levels_SPF0704.dta \
		 

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@


../output/spf_rcons.dta: pull_SPF.py
	python pull_SPF.py

# this is the recipe 
 ../output/levels_SPF0704.dta : graph_SPF_cons_dist.do ../output/spf_rcons.dta $(input)
	stata-se -b run graph_SPF_cons_dist.do
