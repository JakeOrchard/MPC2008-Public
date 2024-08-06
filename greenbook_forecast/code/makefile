## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of processed data files and figures
input = 

output=  ../output/greenbook_rpce.dta \
		../output/greenbook_fore_2008.pdf \
		../output/greenbook08.dta 
		 

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@


../output/greenbook_rpce.dta: pull_greenbook_rpce.py
	python pull_greenbook_rpce.py

# this is the recipe that creates the zip file with cex data  --downloadopt ../input/cexdownloadoptions.yml
../output/greenbook_fore_2008.pdf ../output/greenbook08.dta: graph_greenbook_rpce.do ../output/greenbook_rpce.dta  $(input)
	stata-se -b run graph_greenbook_rpce.do
