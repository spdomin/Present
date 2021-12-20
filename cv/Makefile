#
#  Provides 'make all' to build and 'make clean'
#

all:
	rm -rf *.aux *.nav *.log *.lof *.lot *.out *.snm *.blg *.bbl *.toc *.vrb *.dvi;
	pdflatex --shell-escape dominoCV;

long:
	rm -rf *.aux *.nav *.log *.lof *.lot *.out *.snm *.blg *.bbl *.toc *.vrb *.dvi;
	pdflatex --shell-escape dominoCV_long;

clean:
	-rm -f *.aux *.nav *.log *.lof *.lot *.out *.snm *.blg *.bbl *.toc *.vrb *.dvi;
