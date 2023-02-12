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
Rscript ../graph-cns-tree.R Zm00001eb327910.bestTree Zm00001eb327910_conservedCNSTable.csv Aco003777 Zm00001eb327910 false
