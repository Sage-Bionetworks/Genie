#!/bin/bash
#synapse get syn7537902
#wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
#awk '$3 == "exon" {print}' Homo_sapiens.GRCh37.75.gtf > exon.gtf
#awk '$3 == "gene" {print}' Homo_sapiens.GRCh37.75.gtf > gene.gtf
script_dir=$(dirname $0)

bedtools intersect -a $script_dir/temp.bed  -b $script_dir/exon.gtf -wa | sort | uniq > $script_dir/genie_exons.bed
bedtools intersect -a $script_dir/temp.bed -b $script_dir/exon.gtf -wa -v | sort | uniq > $script_dir/intron_intergenic.bed
bedtools intersect -a $script_dir/temp.bed -b $script_dir/gene.gtf -wa | sort | uniq > $script_dir/gene.bed 
diff $script_dir/gene.bed $script_dir/genie_exons.bed | grep '<' | sed 's/< //' > $script_dir/genie_introns.bed
diff $script_dir/intron_intergenic.bed $script_dir/genie_introns.bed | grep '<' | sed 's/< //'  > $script_dir/genie_intergenic.bed
awk '{print $0"\texon"}' $script_dir/genie_exons.bed > $script_dir/genie_combined.bed
awk '{print $0"\tintron"}' $script_dir/genie_introns.bed >> $script_dir/genie_combined.bed
awk '{print $0"\tintergenic"}' $script_dir/genie_intergenic.bed >> $script_dir/genie_combined.bed
sort -k5 $script_dir/genie_combined.bed > $script_dir/sorted_genie.bed && mv $script_dir/sorted_genie.bed genie_combined.bed
echo -e "Chromosome\tStart_Position\tEnd_Position\tHugo_Symbol\tincludeInPanel\tID\tSEQ_ASSAY_ID\tFeature_Type" > $script_dir/headers.txt
cat $script_dir/headers.txt $script_dir/genie_combined.bed > $script_dir/new_combined.bed
mv $script_dir/new_combined.bed $script_dir/genie_combined.bed
#synapse store genie_combined.bed --id syn7444851