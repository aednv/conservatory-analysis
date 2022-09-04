#!/bin/bash

#BSUB -n 1
#BSUB -J nextflow_job
#BSUB -o nextflow_job.out
#BSUB -e nextflow_job.err
#BSUB -q long
#BSUB -W 30:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"

#To set up R-phytools conda env
#conda config --set channel_priority strict
#conda create --prefix ./phytoolsConda -c conda-forge r-base
#conda activate ./phytoolsConda
#conda install -c conda-forge r-phytools r-codetools r-colorbrewer
#conda deactivate

#set variables
RESUME=false
NO_SEARCH=false
COLORFUL=false

#process input
args=$(getopt --name "$0" --options r:n:c: -- "$@")
eval set -- "$args"

while [[ $# -gt 0 ]]; do
	case "$1" in
		-r) RESUME=$2; shift 2;;
		-n) NO_SEARCH=$2; shift 2;;
		-c) COLORFUL=$2; shift 2;;
		--) shift; break ;;
	esac
done

FASTA=$1
OUTGROUP=$2
REF_GENE=$3

echo "Your input: fasta path-$1, outgroup-$2, refGene-$3. Optional arguments (default false):  resume-$RESUME, noSearch-$NO_SEARCH, colorful-$COLORFUL"  

#only submit job if fasta is found
if [ ! -f $1 ]; then
	echo "Fasta file not found. Example path ./myfasta.fa"
	exit
fi

module load nextflow/20.10.0.5430
if [ $RESUME == false ]; then
	nextflow run cns_tree_generation.nf --startGenes $FASTA --outgroup $OUTGROUP --mainGene $REF_GENE --noSearch $NO_SEARCH --colorful $COLORFUL
fi
if [ $RESUME == true ]; then
	nextflow run cns_tree_generation.nf -resume --startGenes $FASTA --outgroup $OUTGROUP --mainGene $REF_GENE --noSearch $NO_SEARCH --colorful $COLORFUL
fi

