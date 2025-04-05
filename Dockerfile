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


./TrimGalore-0.6.10/trim_galore --fastqc --paired --cores 8 \
./raw/WT_P62/WT_P62_1.fq.gz \
./raw/WT_P62/WT_P62_2.fq.gz \
-o ./trimmed/WT_P62


./TrimGalore-0.6.10/trim_galore --fastqc --paired --cores 8 \
./raw/Y2_P37_28/Y2_P37_28_1.fq.gz \
./raw/Y2_P37_28/Y2_P37_28_2.fq.gz \
-o ./trimmed/Y2_P37_28


./TrimGalore-0.6.10/trim_galore --fastqc --paired --cores 8 \
./raw/YT2_P37_7_17/YT2_P37_7_17_1.fq.gz \
./raw/YT2_P37_7_17/YT2_P37_7_17_2.fq.gz \
-o ./trimmed/YT2_P37_7_17


./STAR/source/STAR --genomeDir ./star_index \
 --runThreadN 20 --readFilesIn \
 ./trimmed/WT_P62/WT_P62_1_val_1.fq.gz \
 ./trimmed/WT_P62/WT_P62_2_val_2.fq.gz \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMattributes Standard \
 --readFilesCommand zcat \
 --outFileNamePrefix ./aligned/WT_P62/


./STAR/source/STAR --genomeDir ./star_index \
 --runThreadN 20 --readFilesIn \
 ./trimmed/Y2_P37_28/Y2_P37_28_1_val_1.fq.gz \
 ./trimmed/Y2_P37_28/Y2_P37_28_2_val_2.fq.gz \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMattributes Standard \
 --readFilesCommand zcat \
 --outFileNamePrefix ./aligned/Y2_P37_28/


./STAR/source/STAR --genomeDir ./star_index \
 --runThreadN 20 --readFilesIn \
 ./trimmed/YT2_P37_7_17/YT2_P37_7_17_1_val_1.fq.gz \
 ./trimmed/YT2_P37_7_17/YT2_P37_7_17_2_val_2.fq.gz \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMattributes Standard \
 --readFilesCommand zcat \
 --outFileNamePrefix ./aligned/YT2_P37_7_17/


featureCounts -T 8 -t exon -g gene_name -s 2 -p --countReadPairs \
-a ./gtf/gencode.v47.annotation.gtf \
-o ./aligned/WT_P62/featureCounts_gene.txt \
./aligned/WT_P62/Aligned.sortedByCoord.out.bam


featureCounts -T 8 -t exon -g gene_name -s 2 -p --countReadPairs \
-a ./gtf/gencode.v47.annotation.gtf \
-o ./aligned/Y2_P37_28/featureCounts_gene.txt \
./aligned/Y2_P37_28/Aligned.sortedByCoord.out.bam


featureCounts -T 8 -t exon -g gene_name -s 2 -p --countReadPairs \
-a ./gtf/gencode.v47.annotation.gtf \
-o ./aligned/YT2_P37_7_17/featureCounts_gene.txt \
./aligned/YT2_P37_7_17/Aligned.sortedByCoord.out.bam


CMD ["/bin/bash"]
