# YAP_TAZ_bulkRNAseq
Bulk RNA sequencing of WT, YAP KO, and YAP/TAZ KO H9 cells

## Getting started

```bash
docker build -t yaptaz .
docker run -it --name ytcontainer yaptaz
```
For running mapping with STAR
# For running in Docker file
# # Get latest STAR source from releases
# RUN git clone https://github.com/alexdobin/STAR.git

# RUN apt-get update && \
#     apt-get install g++ && \
#     apt-get install make && \
#     apt-get install xxd

# # Compile STAR
# RUN cd STAR/source && \
#     make clean && \
#     make STAR

# # Build Genome index.
# RUN cd STAR/source && \
# RUN mkdir -p data/genomes/h38/STAR/ && cd data/genomes/h38/STAR/ && \
#     wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && \
#     wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz && \
#     gunzip *.gz

# #  Make the genome index file. Code from https://github.com/erilu/bulk-rnaseq-analysis
# RUN ./STAR --runThreadN 64 --runMode genomeGenerate --genomeDir ./data/genomes/h38/STAR/ --genomeFastaFiles ./data/genomes/h38/STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ./data/genomes/h38/STAR/Homo_sapiens.GRCh38.99.gtf

To run STAR
```bash
./STAR --help
```
On 64 GB Codespace, takes 40 mins to make the genome index file

Map reads to single sample. Includes unzipping file
```bash
./STAR --runThreadN 64 --genomeDir ./data/genomes/h38/STAR --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn ./WT_P62_1.fq.gz --readFilesCommand zcat --clip3pAdapterSeq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix out_WT_P62_1
```

2-Pass Mapping with STAR
```bash
./STAR --runThreadN 64 \
--runMode alignReads \
--readFilesCommand zcat \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--genomeDir ./data/genomes/h38/STAR \
--readFilesIn ./WT_P62_1.fq.gz ./WT_P62_2.fq.gz \
--outFileNamePrefix out_WT_P62
```

featureCounts Counting reads
Read counts/gene based on mapping results
```bash
featureCounts -s 2 -p -t gene -g gene_id -a ~/dir/annotation.gtf -o counts.txt *.bam
```
Clean the count table
```bash
cut -f1,7-100 counts.txt > featurecounts.txt
```

HISAT2


