# Years for which we download from BLS

all:
	. venv/bin/activate && make -C downloaddata/code
	. venv/bin/activate && make -C download_clean_retail/code
	. venv/bin/activate && make -C appendCEXfiles/code
	. venv/bin/activate && make -C ucccodemappings/code
	. venv/bin/activate && make -C createconsumptionvariables/code
	. venv/bin/activate && make -C checkforsmoothing/code
	. venv/bin/activate && make -C testingconsumptionaggregation/code 
	. venv/bin/activate && make -C processrebatemodule/code
	. venv/bin/activate && make -C createfamilycharacteristics/code
	. venv/bin/activate && make -C nipavariables/code
	. venv/bin/activate && make -C psmjvariables/code
	. venv/bin/activate && make -C psmjsample/code
	. venv/bin/activate && make -C psmjregressions/code
	. venv/bin/activate && make -C descriptive_stats/code
	. venv/bin/activate && make -C graph_nipaVcex/code
	. venv/bin/activate && make -C pcestatistics/code
	. venv/bin/activate && make -C michigan_survey/code
	. venv/bin/activate && make -C survey_prof_fore/code
	. venv/bin/activate && make -C greenbook_forecast/code
	. venv/bin/activate && make -C forecasting/code
	. venv/bin/activate && make -C graph_forecasts_all/code
	. venv/bin/activate && make -C relativepriceofautos/code
	. venv/bin/activate && make -C narrative/code
	. venv/bin/activate && make -C model/code
	. venv/bin/activate && make -C montecarlo/code
	. venv/bin/activate && make -C _finaltablesandfigures/code

# virtual environment
install: venv
	. venv/bin/activate && pip install -r requirements.txt

venv:
	test -d venv || python3 -m venv venv	

# clean build
clean:
	sh cleanbuild.sh
	make all


windows: 
	make -C downloaddata\code
	make -C download_clean_retail\code
	make -C appendCEXfiles\code
	make -C ucccodemappings\code
	make -C createconsumptionvariables\code
	make -C checkforsmoothing\code
	make -C testingconsumptionaggregation/code
	make -C processrebatemodule/code
	make -C createfamilycharacteristics/code
	make -C nipavariables/code
	make -C psmjvariables/code
	make -C psmjsample/code
	make -C psmjregressions/code
	make -C descriptive_stats/code
	make -C graph_nipaVcex/code
	make -C pcestatistics/code
	make -C michigan_survey/code
	make -C forecasting/code
	make -C survey_prof_fore/code
	make -C greenbook_forecast/code
	make -C graph_forecasts_all/code
	make -C relativepriceofautos/code
	make -C narrative/code
	make -C model/code
	make -C montecarlo/code
	make -C _finaltablesandfigures/code
	
	