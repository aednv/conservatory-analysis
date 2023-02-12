#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8192
#SBATCH -p cpu
#SBATCH -t 02:00:00
#SBATCH -o findGenesRoundTwo-%j.out

module load anaconda/2022.10

#load conda environment
conda activate ./../../bioinfToolsConda

#align round 1 Genes
mafft-linsi myRefGeneTest.roundOne.fa > myRefGeneTest.roundOne.aln

#make profile hiddon markov model
hmmbuild -o myRefGeneTest.summary.roundTwo.txt myRefGeneTest.roundTwo.hmm myRefGeneTest.roundOne.aln

#use model to pull top hits out of protein database
hmmsearch -o myRefGeneTest.hmm.roundTwo.txt --noali --tblout myRefGeneTest.search.output.roundTwo.txt myRefGeneTest.roundTwo.hmm ./../../gene_data/protein_database/all.proteins.fasta

#extract all search matches with a higher score than the outgroup
sed '1,3d' myRefGeneTest.search.output.roundTwo.txt > myRefGeneTest.search.output.trimmed.roundTwo.txt

#replace 'Aco003777' with your outgroup gene
cat myRefGeneTest.search.output.trimmed.roundTwo.txt | while read line
        do
                currentGene=$(echo $line | awk '{print $1}')
                echo $currentGene
                samtools faidx ./../../gene_data/protein_database/all.proteins.fasta $currentGene >> myRefGeneTest.roundTwo.fa
                if [ $currentGene == Aco003777 ]; then break; fi
        done
