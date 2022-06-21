# download base ubuntu image (impish)
FROM ubuntu:21.10

# information about this image
LABEL maintainer="adeneve@umass.edu"
LABEL version="0.1"
LABEL description="This is a custom docker image to run Conservatory."

# disable prompt when installing packages
ARG DEBIAN_FRONTEND=noninteractive

# update build tools
RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential
RUN apt-get install -y --no-install-recommends wget

# install blast
RUN apt-get install -y --no-install-recommends ncbi-blast+

# install samtools
RUN apt-get install -y --no-install-recommends samtools

# install bgzip
RUN apt-get install -y --no-install-recommends tabix

# install phyloP (phast)
RUN apt-get install -y --no-install-recommends phast

# install lastz
RUN apt-get install -y --no-install-recommends lastz

# install bedops (wig2bed)
RUN apt-get install -y --no-install-recommends bedops

# install bedtools
RUN apt-get install -y --no-install-recommends bedtools


# install perl modules
RUN apt-get install -y --no-install-recommends libexporter-tiny-perl
RUN apt-get install -y --no-install-recommends liblist-moreutils-perl
RUN apt-get install -y --no-install-recommends libtext-csv-perl
RUN apt-get install -y --no-install-recommends libbio-alignio-stockholm-perl
RUN apt-get install -y --no-install-recommends libdata-dumper-compact-perl
RUN apt-get install -y --no-install-recommends libgetopt-long-descriptive-perl
RUN apt-get install -y --no-install-recommends libstring-diff-perl

# set lang env variable for perl
ENV LANG=C.UTF-8

WORKDIR "/root"
CMD ["/bin/bash"]
