# Use an official Ubuntu base image
FROM ubuntu:latest

# Update and install necessary packages
RUN apt-get update

RUN apt-get install -y perl
RUN apt-get install -y gcc
RUN apt-get install -y g++
RUN apt-get install -y make
RUN apt-get install -y wget
RUN apt-get install -y bzip2
RUN apt-get install -y libcurl4
RUN apt-get install -y libcurl4-openssl-dev

# Installing gtfToGenePred from source
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred \
    && chmod +x gtfToGenePred \
    && mv gtfToGenePred /usr/bin/

# Installing gffread from source
RUN apt-get install -y git
RUN cd /usr/bin \
    && git clone https://github.com/gpertea/gffread.git \
    && cd gffread \
    && make

# Install htslib from source
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libbz2-dev
RUN apt-get install -y liblzma-dev

RUN cd /usr/bin \
    && wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && tar -vxjf htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && make

COPY --chmod=0755 ../scripts/retrieve_seq_from_fasta.pl /usr/local/bin/retrieve_seq_from_fasta.pl
COPY --chmod=0755 ../scripts/convert2annovar.pl /usr/local/bin/convert2annovar.pl
COPY --chmod=0755 ../scripts/coding_change.pl /usr/local/bin/coding_change.pl
COPY --chmod=0755 ../scripts/coding_change_UPDATED.pl /usr/local/bin/coding_change_UPDATED.pl
COPY --chmod=0755 ../scripts/annotate_variation.pl /usr/local/bin/annotate_variation.pl
COPY --chmod=0755 ../scripts/variants_reduction.pl /usr/local/bin/variants_reduction.pl
COPY --chmod=0755 ../scripts/table_annovar.pl /usr/local/bin/table_annovar.pl

# Set environment variables
ENV PATH=/usr/local/bin:/usr/bin/perl:/usr/bin/htslib-1.9:/usr/bin/gtfToGenePred:/usr/bin/gffread:$PATH