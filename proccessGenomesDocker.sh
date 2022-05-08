#!/bin/bash

#BSUB -n 16
#BSUB -R "select[rh=8]"
#BSUB -R "rusage[mem=2048]"
#BSUB -W 16:00
#BSUB -q long
#BSUB -R "span[hosts=1]"
#BSUB -o processGenomesDocker.out
#BSUB -e processGenomesDocker.err

singularity exec --env LANG=C.UTF-8 conservatory-docker.simg perl processGenomes --family Poaceae --verbose --force-orthology --threads 16 --force-run
