## Usage
##  make html -i HIDE=TRUE
###   output: html with hidden content
##  make html -i
###   output: html with revealed content
##  make pdf -i HIDE=TRUE
###   output: pdf with hidden content
##  make pdf -i
###   output: pdf with revealed content

.SUFFIXES: .Rmd .html .rnw .tex .pdf .dvi .ps
HIDE=default
RmdFILES = $(wildcard *.Rmd)
MD = $(patsubst %.Rmd,%.md,$(RmdFILES)) #$(wildcard *.md)
TEX = $(patsubst %.Rmd,%.tex,$(RmdFILES))
HMD = $(patsubst %.Rmd,%.hmd,$(RmdFILES))
DOCX = $(patsubst %.Rmd,%.Dmd,$(RmdFILES))

TARGET = CoralTrends
PANDOC_ARGS=-s -S -i --template='/home/murray/Work/Resources/Scripts/pandoc/template' --bibliography='/home/murray/Work/Resources/References/References.bib' --reference-links 
PANDOC_SC_ARGS =-s -S -i --template='/home/murray/Resources/Scripts/pandoc/template_sc' --bibliography='/home/murray/Work/Resources/References/References.bib'  --reference-links  --self-contained #--csl '/home/murray/Work/Resources/References/oecologia.csl'
PANDOC_XELATEX_ARGS = -s --template='/home/murray/Work/Resources/Scripts/pandoc/default' --bibliography='/home/murray/Work/Resources/References/References.bib' --reference-links -N --pdf-engine=xelatex --toc #--csl '/home/murray/Work/Resources/References/oecologia.csl' 
XELATEX_ARGS=--interaction=nonstopmode --output-driv="xdvipdfmx -vv -V 4"
XELATEX_ARGS=--interaction=batchmode
PANDOC_DOC_ARGS=  -s -S -i --bibliography='/home/murray/Work/Resources/References/References.bib'  #--csl '/home/murray/Work/Resources/References/oecologia.csl' --reference-docx='/home/murray/Work/Resources/Templates/AIMS.docx'

PANDOC=pandoc
PDFS=$(wildcard figures/*.pdf)
JPGS=$(patsubst %.pdf, %.jpg, $(PDFS))
EPSS=$(patsubst %.pdf, %.eps, $(PDFS))
CONVERT = convert -density 300 -resize 33% -background white -flatten

#KNITRMD=$(patsubst %.Rmd, %_knit_.Rmd, $(RmdFILES))
MASTER = CoralTrends.html
DEPENDS = $(TARGET).md

.SECONDARY:  #this is in place to prevent make from removing intermediatory files (like *.Lmd)
%.Dmd: %.Rmd
	echo $@
	$(eval KNITRMD := $(patsubst %.Rmd, %_knit_.dmd, $<))
	cp $< $(KNITRMD)
	sed -i 's/<top>//g' $(KNITRMD)
	sed -i 's/\({r.*\)}/\, dpi=400}/g' $(KNITRMD)
	#echo "library(knitr); purl(\"$<\")" | R --no-save --no-restore
	echo "library(knitr); knit(\"$(KNITRMD)\",output=\"$@\")" | R --no-save --no-restore

docx: $(DOCX)
	$(PANDOC) $(PANDOC_DOC_ARGS) *.Dmd -o $(TARGET).docx

%.Lmd: %.Rmd
	echo $@
	$(eval KNITRMD := $(patsubst %.Rmd, %_knit_.rmd, $<))
	echo "library(knitr); purl(\"$<\")" | R --no-save --no-restore
	cat *.R > all.R
	cp $< $(KNITRMD)
	sed -i "s/='png'/='pdf'/g" $(KNITRMD) #produce pdf graphics
	sed -i "s/='html'/='latex'/g" $(KNITRMD) #produce latex tables
	sed -i "s/sanitize.colnames.function=NULL/sanitize.colnames.function=bold.names/g" $(KNITRMD) #produce latex tables
	sed -i -e '/<top>/{r /home/murray/Work/Resources/Scripts/knitHooksLatex' -e 'd}' $(KNITRMD)
	echo "library(knitr); knit(\"$(KNITRMD)\",output=\"$@\")" | R --no-save --no-restore
	sed -i "s/caption{Table [0-9]*./caption{/g" $@ #remove table counter - latex will do this
	$(eval REPORT_YEAR := $(shell cat parameters/CoralTrends.conf | grep 'final_year' | sed 's/final_year=//g'))
	sed -i "s/final_year/$(REPORT_YEAR)/g" $@ #replace report_year with the value specified in CoralTrends.config
	sed -i "s/coralchange_northern.y1/$(shell cat parameters/CoralTrends.conf | grep 'coralchange_northern.y1' | sed 's/coralchange_northern.y1=//g')/g" $@ #replace report_year with the value specified in CoralTrends.config
	sed -i "s/coralchange_northern.y2/$(shell cat parameters/CoralTrends.conf | grep 'coralchange_northern.y2' | sed 's/coralchange_northern.y2=//g')/g" $@ #replace report_year with the value specified in CoralTrends.config
	sed -i "s/coralchange_central.y1/$(shell cat parameters/CoralTrends.conf | grep 'coralchange_central.y1' | sed 's/coralchange_central.y1=//g')/g" $@ #replace report_year with the value specified in CoralTrends.config
	sed -i "s/coralchange_central.y2/$(shell cat parameters/CoralTrends.conf | grep 'coralchange_central.y2' | sed 's/coralchange_central.y2=//g')/g" $@ #replace report_year with the value specified in CoralTrends.config
	sed -i "s/coralchange_southern.y1/$(shell cat parameters/CoralTrends.conf | grep 'coralchange_southern.y1' | sed 's/coralchange_southern.y1=//g')/g" $@ #replace report_year with the value specified in CoralTrends.config
	sed -i "s/coralchange_southern.y2/$(shell cat parameters/CoralTrends.conf | grep 'coralchange_southern.y2' | sed 's/coralchange_southern.y2=//g')/g" $@ #replace report_year with the value specified in CoralTrends.config

%.tex: %.Lmd
	$(eval TMP := $(patsubst %.Lmd, %.lmd, $<))
	echo $(TMP)
	cp $< $(TMP)
ifeq ($(HIDE),TRUE)
	sed -i '/<div class=\"hidden\".*>/,/^<\/div class=\"hidden\">/d' $(TMP)
endif
	$(PANDOC) $(PANDOC_XELATEX_ARGS) -f 'markdown' *.lmd -o $(TARGET).tex

pdf: $(TEX) #$(TARGET).tex
	sed -i 's/includegraphics{images/includegraphics\[width=\\maxwidth\]{images/g' $(TARGET).tex
	sed -i 's/\\centering/\\centering\\scriptsize/g' $(TARGET).tex #indicate smaller font for tables
	xelatex $(XELATEX_ARGS) $(TARGET).tex
	xelatex $(XELATEX_ARGS) $(TARGET).tex
	xelatex $(XELATEX_ARGS) $(TARGET).tex

%.Hmd: %.Rmd
	$(eval KNITRMD := $(patsubst %.Rmd, %_knit_.rmd, $<))
	echo "library(knitr); purl(\"$<\")" | R --no-save --no-restore
	cat *.R > all.R
	cp $< $(KNITRMD)
	sed -i -e '/<top>/{r /home/murray/Work/Resources/Scripts/knitHooks' -e 'd}' $(KNITRMD)
	echo "library(knitr); knit(\"$(KNITRMD)\",output=\"$@\")" | R --no-save --no-restore

%.hmd: %.Hmd
	$(eval TMP := $(patsubst %.Hmd, %_knit_.hmd, $<))
	cp $< $(TMP)
ifeq ($(HIDE),TRUE)
	sed -i '/<div class=\"hidden\".*>/,/^<\/div class=\"hidden\">/d' $(TMP)
endif

html: $(HMD)
	$(PANDOC) $(PANDOC_ARGS) *.hmd > $(TARGET).html
	sed -i -e '/?config=TeX-AMS_HTML-full/' -e 'd}' *.hmd
	$(PANDOC) $(PANDOC_SC_ARGS) *.hmd > $(TARGET)_sc.html

figs: $(JPGS) $(EPSS)

%.jpg: %.pdf
	@echo ** Building jpg images from pdf versions**
	$(CONVERT) $< $@

%.eps: %.pdf
	pdf2ps -dLanguageLevel=3 $< $@

figures::
	for image_file in $(wildcard figures/*.pdf); \
	do \
	convert -density 300 -quality 100 $${image_file} $${image_file}.jpg; \
	done;

	for image_file in $(wildcard figures/*.pdf); \
	do \
	convert -density 300 -quality 100 $${image_file} $${image_file}.tiff; \
	done;

	for image_file in $(wildcard figures/*.pdf); \
	do \
	pdftops -eps -level3 $${image_file}; \
	done;

clean:
	rm *.toc *.aux *.pdf *.ps *.eps *.log *.lof *.bib *.bbl *.blg *.dvi *.tex *.map *.md
zip:
	zip $(TARGET).zip $(TARGET).html *.R data/*.* figures/*.* scripts/*.* fonts/*.*
