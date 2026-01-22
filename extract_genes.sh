awk -F'\t' '{
    if ($10 ~ /gene_name=/) {
        match($10, /gene_name=([^;]+)/, gname);
        print gname[1];
    }
}' OCC.protein_coding.bed > occ_fusion_genes_extracted.txt
