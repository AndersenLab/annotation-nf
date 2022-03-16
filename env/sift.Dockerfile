FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

RUN conda install -c bioconda sift4g
RUN conda install -c bioconda perl-dbi
RUN conda install -c bioconda perl-bioperl
RUN conda install -c bioconda perl-lwp-simple
RUN apt-get install libswitch-perl

RUN apt-get --allow-releaseinfo-change update && \
   apt-get install -y procps && \
	rm -rf /var/lib/apt/lists/*