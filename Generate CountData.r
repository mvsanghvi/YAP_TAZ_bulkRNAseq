import featureCounts 

#Read counts/gene based on mapping results
featureCounts -s 2 -p -t gene -g gene_id -a ~/dir/annotation.gtf -o counts.txt *.bam
#Clean the count table
cut -f1,7-100 counts.txt > featurecounts.txt