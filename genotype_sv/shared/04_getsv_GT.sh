samples=$(bcftools query -l \
India_China_shared.population_QC.noLR.vcf.gz \
| paste -sd $'\t' -)

printf "CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\t%s\n" "$samples" \
> shared_population_noLR_GT.tsv

bcftools query \
-f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN[\t%GT]\n' \
India_China_shared.population_QC.noLR.vcf.gz \
>> shared_population_noLR_GT.tsv
