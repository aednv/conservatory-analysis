#!/bin/bash

#BSUB -n 1
#BSUB -R "select[rh=8]"
#BSUB -R "rusage[mem=2048]"
#BSUB -W 30:00
#BSUB -q long
#BSUB -R "span[hosts=1]"
#BSUB -o compareCNSDocker.out
#BSUB -e compareCNSDocker.err

singularity exec --env LANG=C.UTF-8 ../conservatory-docker.simg perl compareCNS.pl

