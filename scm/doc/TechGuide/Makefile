#
# Makefile for the CCPP (SCM) Users Guide
#

.PHONY: html latex all clean

all: html latex

latex:
	python -m sphinx -T -b latex -d _build/doctrees -D language=en . latex
	cd latex && latexmk -f -pdf -pdflatex="pdflatex" -use-make index.tex

html:
	python -m sphinx -T -b html -d _build/doctrees -D language=en . html

clean:
	latexmk -CA

