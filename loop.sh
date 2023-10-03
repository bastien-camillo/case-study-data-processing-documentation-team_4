# tmux a
# ctrl + b, d
#!/bin/bash

# Symbolic links
for f in ~/course/data/day2/fastq/PRI-TJPGK*; do ln -s $f .; done

for file in *.fq.gz
do
    vsearch --fastx_uniques $file --fastqout ${file}.vs.fq --minseqlength 30 --strand both
    gzip ${file}.vs.fq
    zcat ${file}.vs.fq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > $file.vs.fq.gz.read_length.txt
    bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U $file.vs.fq.gz --no-unal | samtools view -bS - > $file.bam
done
