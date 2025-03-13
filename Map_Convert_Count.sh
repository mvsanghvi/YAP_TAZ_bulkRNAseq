#Map Reads using HISAT2
# Convert from SAM to BAM

listD = ("WT_P62_1" "WT_P62_2" "Y2_P37_28_1" "Y2_P37_28_2" "YT2_P37_7_17_1" "YT2_P37_7_17_2")
for i in listD:
do
    echo $i
    file_1=$"$i"_1.fq.gz
    file_2=$"$i"_2.fq.gz
    hfile_sam=$"$i"_.sam
    hisat2 -q --dta -p 8 -x grch38/genome -1 ./WT_P62_1.fq.gz -2 ./WT_P62_2.fq.gz  -S WT_P62.sam  -t
    hfile_bam=$"$i"_h.bam
    samtools sort -@ 8 -o hfile_bam hfile_sam

    rm "$hfile_sam"

    yt_folder= $"$i"_all
    stringtie -e -B -p 8 -G grch38/genome ./$hfile_bam -o ./yt_all_gtf/$yt_folder/output_merge.gtf
done