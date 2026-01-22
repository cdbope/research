bcftools view \
  -r chr7 \
  -i 'INFO/SVTYPE="DEL"' \
  gnomad.v4.1.sv.sites.vcf.gz \
  -Oz -o gnomad.chr7.DEL.vcf.gz

