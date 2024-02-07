#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8192
#SBATCH -p cpu
#SBATCH -t 10:00:00
#SBATCH -o buildTree-%j.out

module load anaconda/2022.10

#load conda environment
conda activate ./../../bioinfToolsConda

#remove duplicate proteins
#awk '/^>/ { f = !($0 in a); a[$0]++ } f' myRefGeneTest.roundTwo.fa > myRefGeneTest.roundTwo.fa.trimmed

#align
#mafft-linsi myRefGeneTest.roundTwo.fa.trimmed > myRefGeneTest.roundTwo.fa.trimmed.aln

#build tree
raxml-ng-mpi --all --msa Bradi4g21160.roundTwo.fa.trimmed.aln --model JTT+G --threads 4 --bs-metric fbp,tbe --bs-trees 115
