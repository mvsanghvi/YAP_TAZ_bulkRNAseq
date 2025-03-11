# YAP_TAZ_bulkRNAseq
Bulk RNA sequencing of WT, YAP KO, and YAP/TAZ KO H9 cells

## Getting started

```bash
docker build -t yaptaz .
docker run -it --name ytcontainer yaptaz
```

To run STAR
```bash
./STAR --help
```
On 64 GB Codespace, takes 40 mins to make the genome index file

Map reads to single sample. Includes unzipping file
```bash
./STAR --runThreadN 64 --genomeDir ./data/genomes/h38/STAR --outFilterType BySJout --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn ./WT_P62_1.fq.gz --readFilesCommand zcat --clip3pAdapterSeq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix out_WT_P62_1
```

2-Pass Mapping
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
HTseq Counting Reads
