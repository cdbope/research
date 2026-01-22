install.packages("remotes")
remotes::install_github("FunGeST/Palimpsest")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("PureCN")

PureCN Coverage --out tumor.cov.gz --bam "/media/chbope/Manga/T23-103.merge.bam" --genome hg38 --intervals wgs_10kb.bed
PureCN --out purecn_out --tumor tumor.cov.gz --vcf tumor.vcf.gz --genome hg38
library(data.table)
# pick the table with integer CN per segment (e.g., callAlterations or predictSomatic)
seg <- fread("purecn_out/sample_callAlterations.tsv")
# columns include C (total) and M (minor) → derive nMajor/nMinor
out <- seg[, .(Chromosome=CHR, Start=START, End=END,
               nMajor = pmax(0, C - M), nMinor = M)]

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


Sys.setenv(TMPDIR = "~/mytemp")
dir.create("~/mytemp", showWarnings = FALSE)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ACE")


install.packages("devtools")
devtools::install_github("tgac-vumc/ACE")

install.packages("BiocManager")

R.version.string

Sys.setenv(TMPDIR = "~/mytemp")
dir.create("~/mytemp", showWarnings = FALSE)

###
install.packages(c(
  "caTools",
  "dplyr",
  "ggrepel",
  "tidyr",
  "ggpp"
))

library(ggplot2)
library(caTools)
library(dplyr)
library(ggrepel)
library(tidyr)
library(ggpp)

sample_id <- "T001"
Calls <- read.delim("/home/chbope/extension/nWGS_manuscript_data/data/results/epi2me/epicnv/T001_calls.bed", skip = 1, header = FALSE)
names(Calls) <- c("Chr", "Start", "End", "Bin", "Log2ratio", "Strand")
#print(Calls)

# Read segmentation data
Segs <- read.delim("/home/chbope/extension/nWGS_manuscript_data/data/results/epi2me/epicnv/T001_segs.bed", skip = 1, header = FALSE)
names(Segs) <- c("Chr", "Start", "End", "Bin", "Seg", "Strand")

Calls <- left_join(Calls, Segs)
Calls$Chr <- factor(Calls$Chr, levels = c(1:22, "X", "Y"))

chrom9 <- Calls %>% filter(Chr == 9)

# Read Annot data
Annot <- read.delim("/home/chbope/extension/nWGS_manuscript_data/data/results/cnv/T001_annotatedcnv.csv", header = FALSE)
#print(Annot)
Annot <- Annot %>% 
  separate(V8, c(NA, "Gene"), sep = "=") %>% 
  separate(V6, c(NA, "Score"), sep = "=") %>% 
  separate(V7, c(NA, "LOG2CNT"), sep = "=", convert = TRUE) %>%
  select(c(1,2,6,7,8)) %>%
  filter(Score != "-1") %>%
  filter(Score != "1") %>%
  filter(Gene != "CRLF2")
names(Annot) <- c("Chr", "Start", "Score", "LOG2CNT", "Gene")
Annot$Start <- Annot$Start - 1
Annot$Chr <- as.factor(gsub("chr", "", Annot$Chr))

Annot9 <- Annot %>% filter(Gene == "CDKN2A")

# Merge Calls and Annot
Calls <- left_join(Calls, Annot)

# These chromosomes have very small/absent p-arms that are not represented in the data
p_arm_regions <- data.frame(
  Chr = factor(c("13", "14", "15", "21", "22"), levels = levels(Calls$Chr)),
  Start = c(1, 1, 1, 1, 1),
  End = c(16500000, 16100000, 17083673, 11000000, 13700000),
  Bin = NA,
  Log2ratio = NA,
  Strand = NA,
  Seg = NA,
  Score = NA,
  LOG2CNT = NA,
  Gene = NA
)

# Add the placeholder rows to Calls
Calls <- bind_rows(Calls, p_arm_regions)

# Re-sort by chromosome and position
Calls <- Calls %>% arrange(Chr, Start)

lim <- 3
offset <- 0.1
Calls$Log2ratio_Capped <- ifelse(abs(Calls$Log2ratio) > lim, sign(Calls$Log2ratio)*lim + sign(Calls$Log2ratio)*offset, Calls$Log2ratio)
Calls$Segs_Capped <- ifelse(abs(Calls$Seg) > lim, sign(Calls$Seg)*lim + sign(Calls$Seg)*offset, Calls$Seg)
Calls$flag <- abs(Calls$Log2ratio) > lim
Calls$flagSegs <- abs(Calls$Segs) > lim

# Plot CNV data for all chromosomes
p <- ggplot(Calls) +
  geom_point(data = subset(Calls, !flag), aes(x = Start, y = Log2ratio), colour = "grey", size = 1) +
  geom_point(data = subset(Calls, flag), aes(x = Start, y = Log2ratio_Capped), colour = "red", shape = 17, size = 1) +
  facet_grid(~ Chr, scales = "free_x", space = "free_x") +
  geom_line(data = subset(Calls, !flagSegs), aes(x = Start, y = Seg), colour = "#e41a1c") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_line(colour = "black", size = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_text(size = 12),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "grey", size = 0.4)) +
  scale_y_continuous(minor_breaks = 0)


p2 <- p +
  geom_text_repel(
    aes(x = Start, y = Log2ratio_Capped, label = Gene),
    size = 4,
    color = "black",
    fontface = "bold",
    box.padding = unit(0.5, "lines"),
    #point.padding = unit(0.3, "lines"),
    max.overlaps = 20,       # limit how many labels are shown
    min.segment.length = 0,  # draw connecting lines if needed
    segment.color = "black") + ggtitle(sample_id)
  
  ggsave(filename = out_new_cnv_plot, plot = p2, width = 18, height = 5)
###refit
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("MutationalPatterns")
  