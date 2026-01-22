# Create directory
mkdir -p /humandb/SpliceAI
cd /humandb/SpliceAI

# Download SNV scores + index
wget -c https://storage.googleapis.com/illumina-spliceai/spliceai_scores.masked.snv.hg38.vcf.gz
wget -c https://storage.googleapis.com/illumina-spliceai/spliceai_scores.masked.snv.hg38.vcf.gz.tbi

