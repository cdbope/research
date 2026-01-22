---
output:
  pdf_document:
    fig_height: 3
    latex_engine: pdflatex
    #output: pdf_document
header-includes:
   - \usepackage{helvet}
   - \renewcommand*\familydefault{\sfdefault}
   - \usepackage[T1]{fontenc}
   - \usepackage{xcolor}
   - \usepackage{placeins}
   - \usepackage{float}

geometry: margin=0.5in
---




\begin{flushright}\includegraphics[height=50px]{../../nWGS_manuscript_data/data/log_update} \end{flushright}

# WGS Report - P24 Nanopore Sequencing

## Sample ID: T21-141

##### Read statistics


```
Number of alignments	35727268
% from total alignments	89.69
Number of reads	32106558
Yield [Gb]	205.17
Mean coverage	66.44
Yield [Gb] (>25kb)	0.07
N50	7476
N75	5262
Median length	5352.00
Mean length	6390
Median identity	99.54
Mean identity	109.72
Modal identity	99.6
```


```
The tumor cell content for T21-141 is: 80%
```
### CrossNN Methylation-based Classification

Methylation-based classification is based on **365425** CpG sites (overlap of sites covered in this sample and the model).
At the methylation class (MC) level, the sample has been classified as **glioblastoma, IDH wildtype, subclass RTK II**. 
This prediction has a confidence score of **0.478**.
At the methylation class **family** (MCF) level, the sample has been classified as **glioblastoma, IDH wildtype**.
The MCF prediction has a confidence score of **0.734**.

Scores for the Top 5 entities on MC and MCF level are given below.
Vertical dashed lines indicate the recommended >0.2 cut-off for classification.


\begin{center}\includegraphics[height=150px]{/home/chbope/extension/data/200GMBs/gbm/report/T21-141_markdown_pipeline_report_final_files/figure-latex/unnamed-chunk-4-1} \end{center}


\begin{center}\includegraphics[height=150px]{/home/chbope/extension/data/200GMBs/gbm/report/T21-141_markdown_pipeline_report_final_files/figure-latex/unnamed-chunk-5-1} \end{center}

#### Disclaimer

Methylation-based classification using nanopore whole genome sequencing is a research tool currently under development. It has not been clinically validated in sufficiently large cohorts. Interpretation and implementation of the results in a clinical setting is in the sole responsibility of the treating physician.

### Copy Number Variation Profile

### Full CNV Profile


`Genes annotated in the full CNV profile are amplified (Gain) or deleted (Loss) based on QDNAseq results.`
\includegraphics[height=150px]{../../data/200GMBs/gbm/copy_number_variation/T21-141_cnv_plot_full} 

```
No significant 1p/19q codeletion or Gain 7/Loss 10 detected 
```

```

Details:
* Chr1p: 0/2299 segments show Loss
* Chr19q: 0/576 segments show Loss
* Chr7: 1537/2823 segments show Gain
```

### Copy Number Variation Filter Table

```
The table is filtered for copy number variation events with a score of 2 (amplification) or -2 (homozygous
deletion) and no sex chromosome.
```

\begin{table}[!h]
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}{cccccccc}
\toprule
Chrom & Start & Type & End & SVLEN & Score & LOG2CNT & Gene\\
\midrule
\cellcolor{gray!10}{chr7} & \cellcolor{gray!10}{54900001} & \cellcolor{gray!10}{DUP} & \cellcolor{gray!10}{55600000} & \cellcolor{gray!10}{700000} & \cellcolor{gray!10}{2} & \cellcolor{gray!10}{3.97} & \cellcolor{gray!10}{EGFR}\\
chr12 & 57500001 & DUP & 58050000 & 550000 & 2 & 3.26 & CDK4\\
\bottomrule
\end{tabular}}
\end{table}

### Chromosome 9 CNV/CDK2NA/B Profile


```
Chromosome 9 CNV for visual inspection of CDK2NA/B annotated.
```
\includegraphics[height=150px]{../../data/200GMBs/gbm/copy_number_variation/T21-141_cnv_chr9} 


```
No CDKN2A/B deletions detected
```

### Chromosome 7 CNV/EGFR Profile
```
Chromosome 7 CNV for visual inspection of EGFR annotated.
```
\includegraphics[height=150px]{../../data/200GMBs/gbm/copy_number_variation/T21-141_cnv_chr7} 
```
The EGFR copy number for the sample T21-141 is 39.00
```
```
The vertical red line highlights exons 2-7. Deletion of exons 2-7 results in EGFRviii variant.
```

\begin{center}\includegraphics[width=0.5\linewidth,height=0.3\textheight]{../../data/200GMBs/gbm/tertp/T21-141_egfr_coverage} \end{center}

### IDH1 p.R132 coverage


```
Sequenced coverage of IDH1 p.R132 (red horizontal line) for visual inspection.
```

\begin{center}\includegraphics[width=0.5\linewidth,height=0.3\textheight]{../../data/200GMBs/gbm/tertp/T21-141_idh1_coverage} \end{center}

### TERTp coverage

``` Sequenced coverage of TERTp C228T (left vertical red line) and G250 (right vertical red line) for visual
inspection. ```


\begin{center}\includegraphics[width=0.5\linewidth,height=0.3\textheight]{../../data/200GMBs/gbm/tertp/T21-141_tertp_coverage} \end{center}


### MGMT Methylation Table

```
The table below shows the methylation results for the MGMT gene. Where Mean Methylation Full: All 98 CpG in

the MGMT promoter. Mean Methylation Pyro: Four CpG included in the MGMT pyro kit. Classification by Pyro:

Classification of the pyrosequencing array. Classification by Full: Classification of the full methylation

array.
```

\begin{table}[!h]
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}{ccccc}
\toprule
sample\_id & Mean Methylation Pyro & Mean Methylation Full & Classification by Pyro & Classification by Full\\
\midrule
\cellcolor{gray!10}{T21-141} & \cellcolor{gray!10}{66.0025} & \cellcolor{gray!10}{28.65684} & \cellcolor{gray!10}{Methylated} & \cellcolor{gray!10}{Grey zone, Methylated}\\
\bottomrule
\end{tabular}}
\end{table}
### SNV Calling and Annotation

```
SNVs in the select genes. Only non-synonymous exonic variants that are not known to be benign according to

ClinVar_20240611 are reported. Where GQ: Genotype Quality, Depth: Sequenced depth, AD: Allele Depth, GT:

Allele Genotype, AF: Allele Frequency.
```

\begin{table}[!h]
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}{cccccccccccccccc}
\toprule
Gene.refGene & Chr & Start & End & Ref & Alt & Func.refGene & AAChange.refGene & CLNSIG & COSMIC100 & GQ & Depth & AD & GT & AF & Variant\_caller\\
\midrule
\cellcolor{gray!10}{ID3} & \cellcolor{gray!10}{chr1} & \cellcolor{gray!10}{23559007} & \cellcolor{gray!10}{23559007} & \cellcolor{gray!10}{T} & \cellcolor{gray!10}{C} & \cellcolor{gray!10}{exonic} & \cellcolor{gray!10}{p.T105A} & \cellcolor{gray!10}{0} & \cellcolor{gray!10}{COSV107485314} & \cellcolor{gray!10}{19,28} & \cellcolor{gray!10}{62} & \cellcolor{gray!10}{36,26} & \cellcolor{gray!10}{0/1} & \cellcolor{gray!10}{0.4194} & \cellcolor{gray!10}{Pileup, Merged}\\
PDGFRA & chr4 & 54281602 & 54281602 & C & T & exonic & p.T782M & Conflicting\_classifications\_of\_pathogenicity & COSV57268339 & 20,26 & 71 & 30,40 & 0/1 & 0.5634 & Pileup, Merged\\
\cellcolor{gray!10}{TERT} & \cellcolor{gray!10}{chr5} & \cellcolor{gray!10}{1295113} & \cellcolor{gray!10}{1295113} & \cellcolor{gray!10}{G} & \cellcolor{gray!10}{A} & \cellcolor{gray!10}{upstream} & \cellcolor{gray!10}{No} & \cellcolor{gray!10}{Likely\_pathogenic} & \cellcolor{gray!10}{COSV57196518} & \cellcolor{gray!10}{13,13,34} & \cellcolor{gray!10}{65} & \cellcolor{gray!10}{43,20} & \cellcolor{gray!10}{0/1} & \cellcolor{gray!10}{0.3077} & \cellcolor{gray!10}{Pileup, Merged, ClairS\_TO}\\
KDM4C & chr9 & 7076463 & 7076463 & G & A & exonic & p.S782N & 0 & 0 & 21,110 & 69 & 0,69 & 1/1 & 1.0000 & Pileup, Merged\\
\cellcolor{gray!10}{KDM4C} & \cellcolor{gray!10}{chr9} & \cellcolor{gray!10}{7174673} & \cellcolor{gray!10}{7174673} & \cellcolor{gray!10}{G} & \cellcolor{gray!10}{A} & \cellcolor{gray!10}{exonic} & \cellcolor{gray!10}{p.V602I} & \cellcolor{gray!10}{0} & \cellcolor{gray!10}{COSV67184207} & \cellcolor{gray!10}{20,23} & \cellcolor{gray!10}{60} & \cellcolor{gray!10}{23,36} & \cellcolor{gray!10}{0/1} & \cellcolor{gray!10}{0.6000} & \cellcolor{gray!10}{Pileup, Merged}\\
\addlinespace
TSC1 & chr9 & 132896627 & 132896627 & C & T & exonic & p.G1034S & Conflicting\_classifications\_of\_pathogenicity & COSV53772381 & 20,27 & 69 & 38,31 & 0/1 & 0.4493 & Pileup, Merged\\
\cellcolor{gray!10}{NCOR2} & \cellcolor{gray!10}{chr12} & \cellcolor{gray!10}{124372284} & \cellcolor{gray!10}{124372284} & \cellcolor{gray!10}{C} & \cellcolor{gray!10}{T} & \cellcolor{gray!10}{exonic} & \cellcolor{gray!10}{p.E831K} & \cellcolor{gray!10}{0} & \cellcolor{gray!10}{COSV104422051} & \cellcolor{gray!10}{21,28} & \cellcolor{gray!10}{58} & \cellcolor{gray!10}{27,30} & \cellcolor{gray!10}{0/1} & \cellcolor{gray!10}{0.5172} & \cellcolor{gray!10}{Pileup, Merged}\\
FLT3 & chr13 & 28028284 & 28028284 & T & A & exonic & p.K649N & 0 & 0 & 19,23,31 & 70 & 38,27 & 0/1 & 0.3857 & Pileup, Merged, ClairS\_TO\\
\cellcolor{gray!10}{TP53} & \cellcolor{gray!10}{chr17} & \cellcolor{gray!10}{7673803} & \cellcolor{gray!10}{7673803} & \cellcolor{gray!10}{G} & \cellcolor{gray!10}{A} & \cellcolor{gray!10}{exonic} & \cellcolor{gray!10}{p.R234C} & \cellcolor{gray!10}{Pathogenic/Likely\_pathogenic} & \cellcolor{gray!10}{COSV52662066} & \cellcolor{gray!10}{13,19,34} & \cellcolor{gray!10}{57} & \cellcolor{gray!10}{38,18} & \cellcolor{gray!10}{0/1} & \cellcolor{gray!10}{0.3158} & \cellcolor{gray!10}{Pileup, Merged, ClairS\_TO}\\
CHD7 & chr8 & 60741754 & 60741754 & C & A & exonic & p.P108T & 0 & 0 & 8,23 & 66 & 49,15 & 0/1 & 0.2273 & Merged, ClairS\_TO\\
\addlinespace
\cellcolor{gray!10}{RET} & \cellcolor{gray!10}{chr10} & \cellcolor{gray!10}{43114671} & \cellcolor{gray!10}{43114671} & \cellcolor{gray!10}{G} & \cellcolor{gray!10}{A} & \cellcolor{gray!10}{exonic} & \cellcolor{gray!10}{p.G691S} & \cellcolor{gray!10}{Conflicting\_classifications\_of\_pathogenicity} & \cellcolor{gray!10}{COSV60687096} & \cellcolor{gray!10}{9} & \cellcolor{gray!10}{31} & \cellcolor{gray!10}{4,27} & \cellcolor{gray!10}{0/1} & \cellcolor{gray!10}{0.8710} & \cellcolor{gray!10}{Merged}\\
\bottomrule
\end{tabular}}
\end{table}

\newpage

### Structure Variant Fusion Events

```
No fusion event was detected in this sample.
```

### Disclaimer 

This nanopore whole genome sequencing (nWGS) pipleine is a research tool currently under
development. It has not been clinically validated in sufficiently large cohorts. Interpretation and implementation of the results in a clinical setting is in the sole responsibility of the treating physician.

#### Report generated on 2025-10-08 08:34:05.816012
