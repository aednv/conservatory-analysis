#!/bin/bash

#BSUB -n 16
#BSUB -R "select[rh=8]"
#BSUB -R "rusage[mem=2048]"
#BSUB -W 16:00
#BSUB -q long
#BSUB -R "span[hosts=1]"
#BSUB -o processConservationDocker.out
#BSUB -e processConservationDocker.err

singularity exec --env LANG=C.UTF-8 conservatory-docker.simg perl processConservation --family Poaceae --verbose





