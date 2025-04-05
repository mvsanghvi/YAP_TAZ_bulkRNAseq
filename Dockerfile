FROM python:3


COPY requirements.txt .


RUN pip install -r requirements.txt


RUN mkdir /bulkseq


WORKDIR /bulkseq


# Download human GRCH38 genome index
RUN wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz && \
   gunzip grch38_genome.tar.gz && \
   tar -xvf grch38_genome.tar


# Download human index
RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz


# Download precompiled STAR Index (Human)
RUN wget http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/chrLength.txt \
   http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/chrName.txt \
   http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/chrNameLength.txt \
   http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/chrStart.txt \
   http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/Genome \
   http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/genomeParameters.txt \
   http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/SA \
   http://awspds.refgenie.databio.org/refgenomes.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/star_index__default/SAindex


RUN apt-get update && \
   apt-get install cutadapt fastqc subread


# Install Trim Galore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz && \
   tar xvzf trim_galore.tar.gz

CMD ["/bin/bash"]
