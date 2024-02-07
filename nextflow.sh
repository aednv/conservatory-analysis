#!/bin/bash

#SBATCH -c 1
#SBATCH -p cpu
#SBATCH --mem=8192
#SBATCH -t 24:00:00
#SBATCH -o nextflowJob-%j.out

#To set up phytools / cns mapping conda env
#conda config --set channel_priority strict
#conda create --prefix ./cnsMappingConda -c conda-forge r-base
#conda activate ./cnsMappingConda
#conda install -c conda-forge r-phytools r-codetools r-colorbrewer
#conda deactivate

#To set up bioinformatics tools conda env
#conda create --prefix ./bioinfToolsConda -c bioconda samtools mafft raxml-ng hmmer python pandas

#input argument examples
FASTA=./gene_data/Zm00001eb007950_startGenes.txt
OUTGROUP=Zm00001eb107800_P001
REF_GENE=Zm00001eb007950
PHY_ENV=/work/pi_mbartlett_umass_edu/AmberDeNeve/conservatory-cns-tree-analysis/cnsMappingConda
BIOINF_ENV=/work/pi_mbartlett_umass_edu/AmberDeNeve/conservatory-cns-tree-analysis/bioinfToolsConda

#submit to the cluster
#OUTGROUP=Aco003777 REF_GENE=Zm00001eb007950 FASTA=./gene_data/Zm00001eb007950_startGenes.txt BIOINF_ENV=/work/pi_mbartlett_umass_edu/AmberDeNeve/conservatory-cns-tree-analysis/bioinfToolsConda PHY_ENV=/work/pi_mbartlett_umass_edu/AmberDeNeve/conservatory-cns-tree-analysis/cnsMappingConda <nextflow.sh sbatch

#optional variables (default false)
RESUME=true
NO_SEARCH=false
COLORFUL=false

echo "Your input: fasta path-$FASTA, outgroup-$OUTGROUP, refGene-$REF_GENE. Optional arguments (default false):  resume-$RESUME, noSearch-$NO_SEARCH, colorful-$COLORFUL"  

#only submit job if fasta is found
if [ ! -f $FASTA ]; then
	echo "Fasta file not found. Example path ./myfasta.fa"
	exit
fi

module load nextflow/21.04.3

if [ $RESUME == false ]; then
	nextflow run cns_tree_generation.nf --startGenes $FASTA --outgroup $OUTGROUP --mainGene $REF_GENE --noSearch $NO_SEARCH --colorful $COLORFUL --phytoolsEnv $PHY_ENV --bioinfEnv $BIOINF_ENV
fi
if [ $RESUME == true ]; then
	nextflow run cns_tree_generation.nf -resume --startGenes $FASTA --outgroup $OUTGROUP --mainGene $REF_GENE --noSearch $NO_SEARCH --colorful $COLORFUL --phytoolsEnv $PHY_ENV --bioinfEnv $BIOINF_ENV
fi

