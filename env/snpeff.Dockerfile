# Use an official Ubuntu base image
FROM ubuntu:latest

# Update and install necessary packages
RUN apt-get update

RUN apt-get install -y gcc
RUN apt-get install -y wget
RUN apt-get install -y g++
RUN apt-get install -y make
RUN apt-get install -y bzip2
# install Java Runtime Environment for SnpEff
RUN apt-get install -y default-jre
RUN apt-get install -y unzip


# Existing commands to install htslib dependencies and compile htslib
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev

RUN cd /usr/bin \
    && wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && tar -vxjf htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && make

# Install bcftools - dependent on htslib
RUN cd /usr/bin \
    && wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 \
    && tar -vxjf bcftools-1.9.tar.bz2 \
    && cd bcftools-1.9 \
    && make

# Download and install SnpEff
RUN cd /usr/bin \
    && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
    && unzip snpEff_latest_core.zip -d snpEff \
    && rm snpEff_latest_core.zip

# Set environment variables
ENV PATH=/usr/bin/htslib-1.9:/usr/bin/bcftools-1.9:/usr/bin/snpEff:$PATH