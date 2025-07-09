#
# Makefile for the CCPP (SCM) Users Guide
#

.PHONY: html latex all clean

all: html

latex:
	python -m sphinx -T -b latex -d _build/doctrees -D language=en TechGuide latex
	cd latex && latexmk -f -pdf -pdflatex="pdflatex" -use-make singlecolumnmodelscmforthecommoncommunityphysicspackageccpp.tex

html:
	python -m sphinx -T -b html -d _build/doctrees -D language=en TechGuide html

clean:
	rm -fr _build html latex

cleanlatex:
	cd latex && latexmk -CA
