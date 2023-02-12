#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=8192
#SBATCH -p cpu
#SBATCH -t 02:00:00
#SBATCH -o findGenesRoundOne-%j.out

module load anaconda/2022.10

#load conda environment
conda activate ./../../bioinfToolsConda

#align start Genes. startGenesTest.fa is your fasta with starting protein sequences
mafft-linsi startGenesTest.fa > startGenesTest.aln

#make profile hiddon markov model
hmmbuild -o myRefGeneTest.summary.roundOne.txt myRefGeneTest.roundOne.hmm startGenesTest.aln

#use model to pull top hits out of protein database
hmmsearch -o myRefGeneTest.hmm.roundOne.txt --noali --tblout myRefGeneTest.search.output.roundOne.txt myRefGeneTest.roundOne.hmm ./../../gene_data/protein_database/all.proteins.fasta

#extract all search matches with a higher score than the outgroup
sed '1,3d' myRefGeneTest.search.output.roundOne.txt > myRefGeneTest.search.output.trimmed.roundOne.txt

#replace 'Aco003777' with your outgroup gene
cat myRefGeneTest.search.output.trimmed.roundOne.txt | while read line
	do
		currentGene=$(echo $line | awk '{print $1}')
		echo $currentGene
		samtools faidx ./../../gene_data/protein_database/all.proteins.fasta $currentGene >> myRefGeneTest.roundOne.fa
		if [ $currentGene == Aco003777 ]; then break; fi
	done
