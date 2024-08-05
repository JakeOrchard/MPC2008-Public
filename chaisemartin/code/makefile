## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/psmjsamplemonthly.zip \
		../input/psmjsampleinterview.zip \
		../input/psmjsamplemonthly.dta \
		../input/psmjsampleinterview.dta \
		../input/psmjcleaninterview.dta \
		mpc_utilities.py


# This generates the list of processed data files	
output= \
	..output/test.tex

# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded cex data. if statement checks that the file exists
../input/psmjsampleinterview.zip ../input/psmjsamplemonthly.zip ../input/psmjsampleinterview.dta ../input/psmjsamplemonthly.dta: | ../input
	if [ -e ../../psmjsample/output/$(@F) ] ; then ln -s ../../psmjsample/output/$(@F) $@ ; else exit 1; fi	

../input/psmjcleaninterview.dta: | ../input
	if [ -e ../../psmjregressions/output/$(@F) ] ; then ln -s ../../psmjregressions/output/$(@F) $@ ; else exit 1; fi	

# Create sim link to utility file
mpc_utilities.py:  | ../input
	if [ -e ../../utilities/code/$(@F) ] ; then ln -s ../../utilities/code/$(@F) $@ ; else exit 1; fi	

# this is the recipe that creates the zip file with cex data  --downloadopt ../input/cexdownloadoptions.yml
 ../output/%.tex: %.do $(input)
	stata-se -b run $<
