#Change filename variable for each fastq.gz sample (both forward and reverse)
export filename="YT-P37-7-17-2"

#1. Trim adapters from FASTQ files.
mkdir ./data/trimmed/${filename}
trim_galore --fastqc --paired --cores 8 \
./data/${filename}_L04_R1.fastq.gz \
./data/${filename}_L04_R2.fastq.gz \
-o ./data/trimmed/${filename}

#2. Align data to the reference genome
mkdir ./data/aligned/${filename}/
STAR --genomeDir . \
--runThreadN 20 --readFilesIn \
./data/trimmed/${filename}/${filename}_L04_R1_val_1.fq.gz \
./data/trimmed/${filename}/${filename}_L04_R2_val_2.fq.gz \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--readFilesCommand zcat \
--outFileNamePrefix ./data/aligned/${filename}/

#3. Quantify gene expression
featureCounts -T 8 -t exon -g gene_name -s 2 -p --countReadPairs \
-a ./gencode.v48.annotation.gtf.gz \
-o ./data/aligned/${filename}/featureCounts_exon.txt \
./data/aligned/${filename}/Aligned.sortedByCoord.out.bam