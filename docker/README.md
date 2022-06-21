# Conservatory Docker Image

This environment can be used to run Conservatory, a bioinformatics program for identifying conserved non-coding sequences. Conservatory was developed as a collaboration between the Efroni (HUJI) and Lippman (CSHL) labs. The most recent version of Conservatory can be found here: https://github.com/idanefroni/Conservatory. The docker hub page for this image is [here](https://hub.docker.com/r/aednv/conservatory-docker).

## Included libraries & programs

* Ubuntu 21.10
* NCBI-Blast+
* Samtools
* Phast
* Lastz
* Tabix
* Bedops (wig2bed)
* Bedtools
* Perl modules: List::MoreUtils, Text::CSV, Bio::AlignIO, Data::Dumper, Getopt::Long, String::Diff

## Running the Docker image on a personal computer

Start by installing Docker. Then run this command at your terminal and it will open a bash prompt inside the container.

`docker run -it -v ~/:/host aednv/conservatory-docker`

## Running the Docker image on the MGHPCC cluster (for Bartlett Lab)

The MGHPCC cluster uses Singularity to run docker containers. Detailed info is available [here](https://sylabs.io/guides/2.6/user-guide/singularity_and_docker.html).

* 1. Convert docker image to singularity image. This singularity image may already exist (conservatory-docker.simg), so check before re-building.

      Start interactive session on RedHat 8 node.
      
      `bsub -q interactive -R select[rh=8] -Is /bin/bash`
      
      Build the conservatory-docker.simg image.
      
      `singularity build conservatory-docker.simg docker://aednv/conservatory-docker:latest` 

* 2. Shell into the container while in a RedHat8 interactive session.

      `singularity shell conservatory-docker.simg`
