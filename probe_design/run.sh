# BE SURE TO TYPE "use .vcftools-0.1.14zlib-1.2.6" IF USING BROAD SERVER

#vcf='/broad/mccarroll/dropulation_census_metadata/VCF/CIRMw1w2.clean.corrected.vcf.gz'

# select canonical exons

#awk '{split ($4,a,"_"); {print $1"\t"$2"\t"$3"\t"a[1]"\t"a[3]"\t"$6}}' UCSC_exons.bed > UCSC_exons_modif.bed
#join -1 4 -2 4 <(sort -k4 UCSC_exons_modif.bed ) <(sort -k4 UCSC_canonical.bed) | awk '{print $2"\t"$3"\t"$4"\t"$10"\t"$5"\t"$6}' | bedtools sort -i "-" > UCSC_exons_modif_canonical.bed

# threshold genes

#python threshold_genes.py SW1_SW11_Counts.csv 1000 > high_exp_genes.txt

# create list of positions

#python exons_to_positions.py high_exp_genes.txt UCSC_exons_modif_canonical.bed > high_exp_gene_exons.txt

# get genotypes from VCF 

#vcftools --gzvcf CIRMw1w2w3.no_twins.vcf.gz --out high_exp --012 --positions high_exp_gene_exons.txt --remove-indels &

#python add_gene_exon.py high_exp_gene_exons.txt high_exp.012.pos > high_exp.012.genes

# add probe design info

#python add_mutation.py high_exp.012.genes <(zcat CIRMw1w2w3.no_twins.vcf.gz) > high_exp.012.genes_mutations
cat high_exp.012.genes | awk '{print $1"\t"$2-2"\t"$2+1}' > high_exp.012.snps.bed
bedtools getfasta -fi /mnt/genomes/fasta/hg38/hg38.fa -bed high_exp.012.snps.bed -bedOut > high_exp.012.snps.seq.bed
paste high_exp.012.genes_mutations <(cat high_exp.012.snps.seq.bed | cut -f4) > high_exp.012.genes.full
python add_probe_info.py < high_exp.012.genes.full > high_exp.012.genes.final

bedtools getfasta -fi /mnt/genomes/fasta/hg38/hg38.fa -bed UCSC_exons_modif_canonical.bed -bedOut > UCSC_exons_modif_canonical_seq.bed &
python add_all_regions.py high_exp.012.genes.final2 UCSC_exons_modif_canonical_seq.bed > high_exp.012.genes.regions &
#python add_regions.py final_selection2.txt UCSC_exons_modif_canonical_seq.bed > final_selection_seqs2.txt
python add_melt_and_rep.py < high_exp.012.genes.regions > high_exp.012.final
