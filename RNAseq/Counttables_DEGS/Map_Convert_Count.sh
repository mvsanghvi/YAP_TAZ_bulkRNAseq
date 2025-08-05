#!/bin/bash


# Map Reads using HISAT2
# Convert from SAM to BAM

listD=("WT_P62" "Y2_P37_28" "YT2_P37_7_17")


for i in "${listD[@]}"; do
   echo "$i"
  
   file_1="${i}_1.fq.gz"
   file_2="${i}_2.fq.gz"
   hfile_sam="${i}.sam"
  
   hisat2 -q --dta -p 8 -x grch38/genome -1 "./$file_1" -2 "./$file_2" -S "$hfile_sam" -t
  
   hfile_bam="${i}_h.bam"
  
   samtools sort -@ 8 -o "$hfile_bam" "$hfile_sam"


   rm "$hfile_sam"


   yt_folder=$"$i"_all
   mkdir -p ./yt_all_gtf/$yt_folder
   stringtie -e -B -p 8 -G Homo_sapiens.GRCh38.110.gtf ./$hfile_bam -o ./yt_all_gtf/$yt_folder/output_merge.gtf
done

#Makes gene count matrix
python stringtie-3.0.0.Linux_x86_64/prepDE.py3 -i ./yt_all_gtf -g yt_gene_count_matrix.csv -t yt_transcript_count_matrix.csv
#Makes transcript count matrix: Throwing KeyError in Y2 sample
# /usr/local/python2/bin/python2 stringtie-3.0.0.Linux_x86_64/getTPM.py -i ./yt_all_gtf -g yt_tpm_count_matrix.csv -t yt_tpm_transcript_count_matrix.csv


