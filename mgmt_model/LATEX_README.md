# LaTeX Methods Document - Compilation Guide

## Overview

This directory contains a comprehensive LaTeX version of the MGMT methylation detection pipeline methods documentation, suitable for academic publication in peer-reviewed journals.

## Files

- `METHODS.tex` - Main LaTeX document (774 lines)
- `references.bib` - BibTeX bibliography file (25+ references)
- `Makefile` - Automated compilation script
- `LATEX_README.md` - This guide

## Requirements

### LaTeX Distribution
You need a complete LaTeX installation:

**Ubuntu/Debian:**
```bash
sudo apt-get install texlive-full
```

**CentOS/RHEL:**
```bash
sudo yum install texlive-scheme-full
```

**macOS:**
```bash
brew install mactex
```

**Windows:**
Download and install MiKTeX or TeX Live from their official websites.

### Required LaTeX Packages
The document uses the following packages (included in most full LaTeX distributions):
- `amsmath`, `amsfonts`, `amssymb` - Mathematical typesetting
- `graphicx` - Graphics inclusion
- `natbib` - Bibliography management
- `hyperref` - Hyperlinks and cross-references
- `algorithm`, `algorithmic` - Algorithm typesetting
- `booktabs` - Professional tables
- `geometry` - Page layout
- `listings` - Code listings
- `xcolor` - Color support

## Compilation

### Method 1: Using Makefile (Recommended)

```bash
# Compile complete document with bibliography
make all

# Quick compile without bibliography
make quick

# View the PDF
make view

# Clean auxiliary files
make clean

# Clean all generated files
make cleanall
```

### Method 2: Manual Compilation

```bash
# Full compilation with bibliography
pdflatex METHODS.tex
bibtex METHODS
pdflatex METHODS.tex
pdflatex METHODS.tex

# Quick compilation (no bibliography)
pdflatex METHODS.tex
```

### Method 3: Using LaTeX IDE

Popular LaTeX editors with built-in compilation:
- **TeXstudio** (Cross-platform)
- **Overleaf** (Online)
- **TeXShop** (macOS)
- **WinEdt** (Windows)

## Document Structure

### Main Sections

1. **Abstract** - Comprehensive summary of methods
2. **Introduction** - Background and motivation
3. **Materials and Methods**
   - Data Input and Preprocessing
   - Feature Engineering (98 features)
   - Machine Learning Models (7 algorithms)
   - Model Evaluation
   - Statistical Analysis
   - Computational Implementation
4. **Software Implementation**
5. **Quality Control and Limitations**
6. **Reproducibility**
7. **Conclusion**
8. **References** (25+ peer-reviewed citations)

### Mathematical Notation

The document includes comprehensive mathematical formulations:

- **Feature Engineering**: 30+ mathematical equations
- **Machine Learning**: Algorithm-specific formulations
- **Statistical Analysis**: Hypothesis testing, FDR correction
- **Evaluation Metrics**: ROC-AUC, precision, recall formulations

### Key Features

- **Publication-ready formatting**
- **Comprehensive mathematical notation**
- **Algorithm pseudocode**
- **Professional bibliography**
- **Cross-references and hyperlinks**
- **Standardized academic structure**

## Customization

### Journal-Specific Formatting

The document can be easily adapted for different journals:

**Nature/Science format:**
```latex
\documentclass[10pt,twocolumn]{article}
\usepackage[margin=2cm]{geometry}
```

**IEEE format:**
```latex
\documentclass[conference]{IEEEtran}
```

**Springer format:**
```latex
\documentclass{svjour3}
```

### Bibliography Styles

Change bibliography style in the main document:
```latex
\bibliographystyle{nature}     % Nature journals
\bibliographystyle{ieee}       % IEEE journals
\bibliographystyle{springer}   % Springer journals
\bibliographystyle{apa}        % APA style
```

## Output

The compilation produces:
- `METHODS.pdf` - Final formatted document
- `METHODS.aux` - Auxiliary file for cross-references
- `METHODS.bbl` - Processed bibliography
- `METHODS.log` - Compilation log

## Troubleshooting

### Common Issues

**1. Missing Packages**
```
! LaTeX Error: File 'package.sty' not found
```
Solution: Install missing package or use full LaTeX distribution.

**2. Bibliography Issues**
```
LaTeX Warning: Citation 'key' on page X undefined
```
Solution: Run `bibtex` step in compilation sequence.

**3. Mathematical Symbols**
```
! Undefined control sequence
```
Solution: Ensure `amsmath` and `amsfonts` packages are loaded.

### Compilation Log Analysis

Check `METHODS.log` for:
- Warning messages
- Missing references
- Overfull/underfull boxes
- Package loading issues

### File Permissions

Ensure write permissions for LaTeX auxiliary files:
```bash
chmod 755 *.tex *.bib
```

## Version Control

For collaborative editing:
```bash
# Track only source files
git add METHODS.tex references.bib Makefile
git add LATEX_README.md

# Ignore generated files
echo "*.aux" >> .gitignore
echo "*.bbl" >> .gitignore
echo "*.blg" >> .gitignore
echo "*.log" >> .gitignore
echo "*.pdf" >> .gitignore
```

## Integration with Manuscript Preparation

### Embedding in Larger Documents

The methods section can be included in a larger manuscript:
```latex
\input{METHODS}
```

### Converting to Other Formats

**Word Document:**
```bash
pandoc METHODS.tex -o METHODS.docx
```

**HTML:**
```bash
pandoc METHODS.tex -o METHODS.html --mathjax
```

## Quality Assurance

### Pre-submission Checklist

- [ ] All equations numbered and referenced
- [ ] All citations properly formatted
- [ ] Mathematical notation consistent
- [ ] Algorithms clearly presented
- [ ] No compilation warnings
- [ ] Bibliography complete and accurate
- [ ] Cross-references working
- [ ] Page layout appropriate for target journal

### Validation

Run these commands to validate document quality:
```bash
# Check for undefined references
grep -n "LaTeX Warning" METHODS.log

# Validate bibliography
bibtex METHODS 2>&1 | grep -i error

# Check mathematical equations
grep -n "equation" METHODS.tex
```

## Academic Standards

The document meets requirements for:
- **High-impact journals** (Nature, Science, Cell)
- **Specialized journals** (Bioinformatics, BMC, PLOS)
- **Conference proceedings** (RECOMB, ISMB, ICML)
- **Thesis chapters**
- **Grant applications**

## Support

For LaTeX-specific issues:
- [LaTeX Stack Exchange](https://tex.stackexchange.com/)
- [Overleaf Documentation](https://www.overleaf.com/learn)
- [CTAN Package Repository](https://ctan.org/)

For document content:
- Refer to the main pipeline documentation
- Check the source code comments
- Review the comprehensive README files

---

**Note**: This LaTeX document provides the mathematical rigor and professional formatting required for peer-reviewed publication of computational methods in bioinformatics and machine learning.

