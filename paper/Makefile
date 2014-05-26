lighter_paper.pdf: lighter_paper.tex lighter_paper.bib
	latex lighter_paper
	bibtex lighter_paper
	latex lighter_paper
	bibtex lighter_paper
	latex lighter_paper
	#dvipdf lighter_paper.dvi

.PHONY: clean
clean:
	rm -rf *.bbl *.blg *.log *.out *.pdf *.aux
