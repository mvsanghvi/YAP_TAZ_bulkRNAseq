FROM python:3

# Install conda with libraries for RNA-seq analysis
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
RUN bash Miniforge3-$(uname)-$(uname -m).sh
RUN conda create -y -n rnaseq_env python=3.9 && \
   conda activate rnaseq_env && \
   conda install -y -c bioconda -c conda-forge \
       fastqc \
       fastq-screen \
       multiqc \
       trim-galore \
       cutadapt \
       star=2.7.1a \
       bowtie2 \
       samtools \
       subread \
       sra-tools


RUN mkdir /bulkseq


WORKDIR /bulkseq


# Download human index: https://www.gencodegenes.org/human/
RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz




# Download precompiled STAR Index (Human)
RUN wget http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/chrLength.txt \
  http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/chrName.txt \
  http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/chrNameLength.txt \
  http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/chrStart.txt \
  http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/Genome \
  http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/genomeParameters.txt \
  http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/SA \
  http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/SAindex


CMD ["/bin/bash"]