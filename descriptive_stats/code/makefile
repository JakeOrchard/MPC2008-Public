## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/psmjsampleinterview.dta 


# This generates the list of processed data files	
output= \
		../output/interview_randomnessPB.tex \
		../output/interview_randomness.tex

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded cex data. if statement checks that the file exists
../input/psmjsampleinterview.dta: | ../input
	if [ -e ../../psmjsample/output/$(@F) ] ; then ln -s ../../psmjsample/output/$(@F) $@ ; else exit 1; fi	



# we first delete the output folders to make sure all data is new 
../output/interview_randomnessPB.tex ../output/interview_randomness.tex: schedule_table.do  $(input)
	stata-se -b run schedule_table.do
	rm -f *.log
