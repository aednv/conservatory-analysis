#!/bin/bash

#BSUB -n 1
#BSUB -R "select[rh=8]"
#BSUB -R "rusage[mem=2048]"
#BSUB -W 30:00
#BSUB -q long
#BSUB -R "span[hosts=1]"
#BSUB -o makeCNSgffDocker.out
#BSUB -e makeCNSgffDocker.err

singularity exec --env LANG=C.UTF-8 ../conservatory-docker.simg perl -w scripts/buildConCNS --family Poaceae --format GFF --CNSfile CNS/Poaceae.csv > Poaceae.cns.gff






