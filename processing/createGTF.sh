#!/bin/bash
script_dir=$(dirname $0)
cd $script_dir
wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip $script_dir/Homo_sapiens.GRCh37.75.gtf.gz
awk '$3 == "exon" {print}' $script_dir/Homo_sapiens.GRCh37.75.gtf > $script_dir/exon.gtf
awk '$3 == "gene" {print}' $script_dir/Homo_sapiens.GRCh37.75.gtf > $script_dir/gene.gtf



