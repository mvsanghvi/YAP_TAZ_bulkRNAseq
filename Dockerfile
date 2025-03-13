FROM python:3

COPY requirements.txt .

RUN pip install -r requirements.txt

RUN mkdir /bulkseq

WORKDIR /bulkseq

#Install HISAT [HISAT2](https://daehwankimlab.github.io/hisat2/manual/)
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && \
    unzip download
    # cp ./hisat2-2.2.1/* /usr/bin/

# Download human GRCH38 genome index

RUN wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz && \
    gunzip grch38_genome.tar.gz && \
    tar -xvf grch38_genome.tar

# Install samtools: https://github.com/samtools/samtools/blob/develop/INSTALL

RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2


RUN bunzip2 samtools-1.21.tar.bz2 && \
    tar -xvf samtools-1.21.tar && \
    cd samtools-1.21 && \
    ./configure && \
    make && \
    make install

# Install stringtie: https://github.com/gpertea/stringtie
RUN wget https://ccb.jhu.edu/software/stringtie/dl/stringtie-3.0.0.Linux_x86_64.tar.gz && \
    gunzip stringtie-3.0.0.Linux_x86_64.tar.gz && \
    tar -xvf stringtie-3.0.0.Linux_x86_64.tar
CMD ["/bin/bash"]
