#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=8192
#SBATCH -p cpu
#SBATCH -t 02:00:00
#SBATCH -o mapCnsData-%j.out

module load anaconda/2022.10

#load conda environment
conda activate ./../../bioinfToolsConda

#extract just gene names from list
cat Bradi4g21160.roundTwo.fa.trimmed.aln | grep '>' | cut -b 2- > Bradi4g21160_TreeGenes.txt

#replace with your gene id
python ../csv-database-cns-tree-generation.py Bradi4g21160 ../../gene_data
