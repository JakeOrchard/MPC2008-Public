## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of processed data files and figures
input = 

output=    ../output/all_ffA.eps

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@


../output/all_ffA.eps:  montecarlo_lagexpenditure.do
	stata-se -b run montecarlo_lagexpenditure.do
