FROM python:3

# Get latest STAR source from releases
RUN git clone https://github.com/alexdobin/STAR.git

RUN apt-get update && \
    apt-get install g++ && \
    apt-get install make && \
    apt-get install xxd

# Compile STAR
RUN cd STAR/source && \
    make clean && \
    make STAR

# Build Genome index
RUN wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
    wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz && \
    gunzip *.gz
RUN mkdir -p data/genomes/h38/STAR/
RUN mv ./Homo_sapiens.GRCh38.99.gtf ./data/genomes/h38/STAR/
RUN mv ./Homo_sapiens.GRCh38.dna.primary_assembly.fa ./data/genomes/h38/STAR/
# RUN ./STAR --runThreadN 64 --runMode genomeGenerate --genomeDir ./data/genomes/h38/STAR/ --genomeFastaFiles ./data/genomes/h38/STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ./data/genomes/h38/STAR/Homo_sapiens.GRCh38.99.gtf

WORKDIR /STAR/source

CMD ["/bin/bash"]
