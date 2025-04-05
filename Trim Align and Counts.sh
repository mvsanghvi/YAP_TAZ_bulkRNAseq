./TrimGalore-0.6.10/trim_galore --fastqc --paired --cores 8 \
./raw/WT_P62/WT_P62_1.fq.gz \
./raw/WT_P62/WT_P62_2.fq.gz \
-o ./trimmed/WT_P62

./STAR/source/STAR --genomeDir ./star_index \
 --runThreadN 20 --readFilesIn \
 ./trimmed/WT_P62/WT_P62_1_val_1.fq.gz \
 ./trimmed/WT_P62/WT_P62_2_val_2.fq.gz \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMattributes Standard \
 --readFilesCommand zcat \
 --outFileNamePrefix ./aligned/WT_P62/

featureCounts -T 8 -t exon -g gene_name -s 2 -p --countReadPairs \
-a ./gtf/gencode.v47.annotation.gtf \
-o ./aligned/WT_P62/featureCounts_gene.txt \
./aligned/WT_P62/Aligned.sortedByCoord.out.bam

#Download featureCounts and delete everything else to make room for next batch

./TrimGalore-0.6.10/trim_galore --fastqc --paired --cores 8 \
./raw/Y2_P37_28/Y2_P37_28_1.fq.gz \
./raw/Y2_P37_28/Y2_P37_28_2.fq.gz \
-o ./trimmed/Y2_P37_28

./STAR/source/STAR --genomeDir ./star_index \
 --runThreadN 20 --readFilesIn \
 ./trimmed/Y2_P37_28/Y2_P37_28_1_val_1.fq.gz \
 ./trimmed/Y2_P37_28/Y2_P37_28_2_val_2.fq.gz \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMattributes Standard \
 --readFilesCommand zcat \
 --outFileNamePrefix ./aligned/Y2_P37_28/

featureCounts -T 8 -t exon -g gene_name -s 2 -p --countReadPairs \
-a ./gtf/gencode.v47.annotation.gtf \
-o ./aligned/Y2_P37_28/featureCounts_gene.txt \
./aligned/Y2_P37_28/Aligned.sortedByCoord.out.bam

#Download featureCounts and delete everything else to make room for next batch

./TrimGalore-0.6.10/trim_galore --fastqc --paired --cores 8 \
./raw/YT2_P37_7_17/YT2_P37_7_17_1.fq.gz \
./raw/YT2_P37_7_17/YT2_P37_7_17_2.fq.gz \
-o ./trimmed/YT2_P37_7_17


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
-o ./aligned/YT2_P37_7_17/featureCounts_gene.txt \
./aligned/YT2_P37_7_17/Aligned.sortedByCoord.out.bam

#Download featureCounts and delete everything else to make room for next batch