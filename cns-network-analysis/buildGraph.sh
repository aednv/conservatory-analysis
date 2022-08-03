#!/bin/bash

#BSUB -n 2
#BSUB -R "select[rh=8]"
#BSUB -R "rusage[mem=20000]"
#BSUB -W 40:00
#BSUB -q long
#BSUB -R "span[hosts=1]"
#BSUB -o buildGraph.out
#BSUB -e buildGraph.err

singularity exec ../../singularity-container-images/networkx-docker.simg python buildGraph.py

