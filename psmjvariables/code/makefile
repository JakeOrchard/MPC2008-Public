## Definitions

# Folders we will be looking for in the task directory
folders := ../input ../output

# This generates the list of zip files we should see in the input folder 	
input= \
		../input/fmliquarterly.parquet \
		../input/mtbimonthly.parquet \
		../input/cohortvariables.parquet \
		../input/interviewvariables.parquet \
		../input/rebatecuid.parquet \
		../input/rebateinterview.parquet \
		../input/rebatemonthly.parquet \
		../input/nipavariablesinterview.parquet \
		../input/nipavariablesmonthly.parquet \


# This generates the list of processed data files	
output= \
		../output/psmjexpendituresmonthly.parquet \
		../output/psmjexpendituresinterview.parquet \
		../output/psmjvariablesinterview.parquet \
		../output/psmjvariablesmonthly.parquet 


# Recipes
all: $(folders) $(input) $(output)

# Create folders if they do not already exist
$(folders):
	mkdir $@

# Create sim link to file with downloaded cex data. if statement checks that the file exists
../input/fmliquarterly.parquet ../input/mtbimonthly.parquet: | ../input
	if [ -e ../../createconsumptionvariables/output/$(@F) ] ; then ln -s ../../createconsumptionvariables/output/$(@F) $@ ; else exit 1; fi	

../input/cohortvariables.parquet  ../input/interviewvariables.parquet:  | ../input
	if [ -e ../../createfamilycharacteristics/output/$(@F) ] ; then ln -s ../../createfamilycharacteristics/output/$(@F) $@ ; else exit 1; fi	

../input/rebatecuid.parquet  ../input/rebateinterview.parquet ../input/rebatemonthly.parquet:  | ../input
	if [ -e ../../processrebatemodule/output/$(@F) ] ; then ln -s ../../processrebatemodule/output/$(@F) $@ ; else exit 1; fi	

../input/nipavariablesinterview.parquet ../input/nipavariablesmonthly.parquet :  | ../input
	if [ -e ../../nipavariables/output/$(@F) ] ; then ln -s ../../nipavariables/output/$(@F) $@ ; else exit 1; fi	


# this is the recipe
../output/psmjexpendituresinterview.parquet ../output/psmjexpendituresmonthly.parquet: psmjexpenditures.py ../input/fmliquarterly.parquet ../input/mtbimonthly.parquet
	rm -r ../output
	mkdir ../output
	python $<

../output/psmjvariablesinterview.parquet ../output/psmjvariablesmonthly.parquet: psmjmergedata.py ../output/psmjexpendituresinterview.parquet ../output/psmjexpendituresmonthly.parquet $(input)
	python $<
