# Folders for storing tables and figures
folders := tables figures

# PSMJ regressions
REGTABLES := $(wildcard ../psmjregressions/output/*.tex)
REGTABLELINKS := $(patsubst ../psmjregressions/output/%.tex, tables/%.tex, $(REGTABLES))

# SA figures
SAFIGURES := $(wildcard ../sunabraham/output/*.eps)
SAFIGURELINKS := $(patsubst ../sunabraham/output/%.eps, figures/%.eps, $(SAFIGURES))

SATABLES := $(wildcard ../sunabraham/output/*.tex)
SATABLELINKS := $(patsubst ../sunabraham/output/%.tex, tables/%.tex, $(SATABLES))

#Descriptive Stats
DESTABLES :=  $(wildcard ../descriptive_stats/output/*.tex)
DESTABLELINKS := $(patsubst ../descriptive_stats/output/%.tex, tables/%.tex, $(DESTABLES))

#PCE v CEX
PCETABLES :=  $(wildcard ../graph_nipaVcex/output/*.tex)
PCETABLELINKS := $(patsubst ../graph_nipaVcex/output/%.tex, tables/%.tex, $(PCETABLES))

PCEFIGURES  :=  $(wildcard ../graph_nipaVcex/output/*.pdf)
PCEFIGURELINKS  := $(patsubst ../graph_nipaVcex/output/%.pdf, figures/%.pdf, $(PCEFIGURES))


all: logbook.pdf $(folders) $(REGTABLELINKS) $(SAFIGURELINKS) $(SATABLELINKS) $(DESTABLELINKS) 

$(folders):
	mkdir $@ 

$(REGTABLELINKS): 
	ln -s ../../psmjregressions/output/$(@F) $@

$(DESTABLELINKS): 
	ln -s ../../descriptive_stats/output/$(@F) $@

$(SAFIGURELINKS) $(SATABLELINKS): 
	ln -s ../../sunabraham/output/$(@F) $@	

$(PCEFIGURELINKS) $(PCETABLELINKS): 
	ln -s ../../graph_nipaVcex/output/$(@F) $@	

logbookentries = $(shell grep '\\input{.*entries' logbook.tex | sed 's/^.*\\input{\(.*\.tex\).*/\1/')

logbook.pdf: logbook.tex $(logbookentries) $(REGTABLELINKS) $(SAFIGURELINKS) $(SATABLELINKS) $(DESTABLELINKS) $(PCEFIGURELINKS) $(PCETABLELINKS)
	$(shell cat tables/scalar*.tex > tables/all_scalars.tex)
	pdflatex -shell-escape -draftmode logbook.tex
	pdflatex -shell-escape logbook.tex
	rm logbook.log logbook.out logbook.toc logbook.aux
