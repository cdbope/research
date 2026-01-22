# Makefile for MGMT Methods LaTeX Document

# Main document
MAIN = METHODS
TEX_FILE = $(MAIN).tex
PDF_FILE = $(MAIN).pdf
BIB_FILE = references.bib

# LaTeX compiler
LATEX = pdflatex
BIBTEX = bibtex

# Default target
all: $(PDF_FILE)

# Compile PDF
$(PDF_FILE): $(TEX_FILE) $(BIB_FILE)
	$(LATEX) $(TEX_FILE)
	$(BIBTEX) $(MAIN)
	$(LATEX) $(TEX_FILE)
	$(LATEX) $(TEX_FILE)

# Quick compile (no bibliography)
quick: $(TEX_FILE)
	$(LATEX) $(TEX_FILE)

# Clean auxiliary files
clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz

# Clean all generated files
cleanall: clean
	rm -f $(PDF_FILE)

# View PDF (if available)
view: $(PDF_FILE)
	if command -v evince > /dev/null; then \
		evince $(PDF_FILE) & \
	elif command -v okular > /dev/null; then \
		okular $(PDF_FILE) & \
	elif command -v xpdf > /dev/null; then \
		xpdf $(PDF_FILE) & \
	else \
		echo "No PDF viewer found. Please install evince, okular, or xpdf."; \
	fi

# Help
help:
	@echo "Available targets:"
	@echo "  all      - Compile the complete document with bibliography"
	@echo "  quick    - Quick compile without bibliography"
	@echo "  clean    - Remove auxiliary files"
	@echo "  cleanall - Remove all generated files including PDF"
	@echo "  view     - Open the PDF in a viewer"
	@echo "  help     - Show this help message"

# Force targets
.PHONY: all quick clean cleanall view help

