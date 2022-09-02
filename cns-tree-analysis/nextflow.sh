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

module load nextflow/20.10.0.5430

nextflow run cns_tree_generation.nf -resume
