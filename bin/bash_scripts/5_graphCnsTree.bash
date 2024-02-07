#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=8192
#SBATCH -p cpu
#SBATCH -t 01:00:00
#SBATCH -o graphCnsTree-%j.out

module load anaconda/2022.10

#load conda environment
conda activate ./../../cnsMappingConda

#replace with your tree file, cns mapping file, outgroup gene, gene_id, and color preference
Rscript ../graph-cns-tree.R Bradi4g21160.roundTwo.fa.trimmed.aln.raxml.supportFBP Bradi4g21160_conservedCNSTable.csv Bradi4g40270 Bradi4g21160 false
