lighter_paper.pdf: lighter_paper.tex lighter_paper.bib
	pdflatex lighter_paper
	bibtex lighter_paper
	pdflatex lighter_paper
	bibtex lighter_paper
	pdflatex lighter_paper

.PHONY: clean
clean:
	rm -rf *.bbl *.blg *.log *.out *.pdf *.aux
